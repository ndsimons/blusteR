#' Discover enriched CDR3 motifs (local similarity)
#'
#' Analogous to the GLIPH2 local-motif step.  Extracts k-mers from CDR3H
#' (and optionally CDR3L) sequences, then tests each k-mer for
#' enrichment against the OAS background using a one-sided binomial test.
#'
#' Because BCR CDR3H regions are longer and more diverse than TCR CDR3β
#' (median ~15 aa vs ~13 aa), blusteR uses k = 4 and k = 5 by default
#' (GLIPH2 uses k = 2–3 for "local convergence" motifs).  Somatic
#' hypermutation is handled by also testing degenerate motifs where one
#' position is replaced with a physico-chemical group wildcard.
#'
#' @param bcr_data A \code{data.table} in blusteR standard format.
#' @param reference A \code{bluster_reference} object, or \code{NULL} to
#'   auto-build one.
#' @param chain Which chain to analyse: \code{"heavy"} (default),
#'   \code{"light"}, or \code{"both"}.
#' @param k_sizes Integer vector of k-mer lengths (default \code{c(4, 5)}).
#' @param min_freq Minimum number of input sequences that must contain a
#'   motif for it to be tested (default 3).
#' @param max_freq Maximum fraction of input sequences a motif may appear
#'   in (default 0.20).  Motifs found in more than this fraction are
#'   germline/framework-common rather than specificity signals and are
#'   dropped to avoid merging unrelated cells into one giant cluster.
#' @param p_cutoff Adjusted p-value threshold (default 0.05).
#' @param shm_degenerate Logical; also test degenerate motifs with one
#'   position masked to its physico-chemical group (default TRUE).
#' @param p_adjust_method Method for \code{p.adjust} (default \code{"BH"}).
#'
#' @return A \code{data.table} with columns: \code{motif, k, chain,
#'   count_obs, count_ref, fold_enrichment, pvalue, pvalue_adj, member_ids}.
#'
#' @export
bluster_motifs <- function(bcr_data,
                         reference = NULL,
                         chain = c("heavy", "light", "both"),
                         k_sizes = .DEFAULT_KMER_SIZES,
                         min_freq = .DEFAULT_LOCAL_MIN_FREQ,
                         max_freq = .DEFAULT_LOCAL_MAX_FREQ,
                         p_cutoff = .DEFAULT_PVALUE,
                         shm_degenerate = TRUE,
                         p_adjust_method = "BH") {

  chain <- match.arg(chain)

  # Build reference if not provided
  if (is.null(reference)) {
    message("[blusteR] No reference provided; building default reference...")
    reference <- build_reference()
  }

  # Select chains to analyse
  chains_to_test <- switch(chain,
    heavy = "heavy",
    light = "light",
    both  = c("heavy", "light")
  )

  all_results <- list()

  for (ch in chains_to_test) {
    seqs <- bcr_data[chain_type == ch, .(
      seq_id  = cell_id,
      cdr3    = cdr3_aa,
      v_gene  = v_gene
    )]
    seqs <- seqs[!is.na(cdr3) & nchar(cdr3) >= max(k_sizes)]

    if (nrow(seqs) == 0) {
      message("[blusteR] No ", ch, " chain sequences to test.")
      next
    }

    message(sprintf("[blusteR] Testing local motifs in %d %s chain CDR3s...",
                    nrow(seqs), ch))

    for (k in k_sizes) {
      # Extract k-mers from input sequences
      kmer_map <- .extract_kmers_per_seq(seqs$cdr3, k)
      kmer_counts <- table(unlist(kmer_map))

      # Filter to minimum frequency, and drop promiscuous (too-common) motifs
      max_count <- max_freq * nrow(seqs)
      candidates <- names(kmer_counts)[kmer_counts >= min_freq &
                                       kmer_counts <= max_count]

      if (length(candidates) == 0) next

      # Get background frequencies
      bg_key <- as.character(k)
      bg <- reference$kmer_freq[[bg_key]]

      if (is.null(bg)) {
        message("[blusteR] No background k-mer frequencies for k=", k,
                "; skipping.")
        next
      }

      # Test each candidate motif
      results <- .test_motif_enrichment(
        candidates  = candidates,
        kmer_map    = kmer_map,
        seq_ids     = seqs$seq_id,
        bg          = bg,
        n_input     = nrow(seqs),
        k           = k,
        chain_label = ch
      )

      all_results[[paste0(ch, "_k", k)]] <- results

      # Optionally test SHM-degenerate motifs
      if (shm_degenerate) {
        degen_motifs <- .generate_degenerate_motifs(candidates, k)
        if (length(degen_motifs) > 0) {
          degen_map <- .extract_degenerate_kmers_per_seq(seqs$cdr3, degen_motifs)
          degen_counts <- .count_degenerate_motifs(degen_motifs, degen_map)
          degen_cands <- names(degen_counts)[degen_counts >= min_freq &
                                             degen_counts <= max_count]
          if (length(degen_cands) > 0) {
            degen_results <- .test_degenerate_enrichment(
              candidates  = degen_cands,
              degen_map   = degen_map,
              seq_ids     = seqs$seq_id,
              bg          = bg,
              n_input     = nrow(seqs),
              k           = k,
              chain_label = ch
            )
            all_results[[paste0(ch, "_k", k, "_degen")]] <- degen_results
          }
        }
      }
    }
  }

  if (length(all_results) == 0) {
    message("[blusteR] No enriched motifs found.")
    return(data.table::data.table(
      motif = character(0), k = integer(0), chain = character(0),
      count_obs = integer(0), freq_obs = numeric(0),
      freq_ref = numeric(0), fold_enrichment = numeric(0),
      pvalue = numeric(0), pvalue_adj = numeric(0),
      member_ids = list()
    ))
  }

  out <- data.table::rbindlist(all_results, fill = TRUE)

  # Multiple-testing correction across all motifs
  out[, pvalue_adj := p.adjust(pvalue, method = p_adjust_method)]

  # Filter to significant
  out <- out[pvalue_adj < p_cutoff]

  # Sort by enrichment
  data.table::setorder(out, pvalue_adj)

  message(sprintf("[blusteR] Found %d enriched motifs.", nrow(out)))
  out[]
}


# ---- internal motif helpers ------------------------------------------

#' Extract all k-mers from each sequence
#' @keywords internal
.extract_kmers_per_seq <- function(cdr3_vec, k) {

  lapply(cdr3_vec, function(s) {
    n <- nchar(s)
    if (n < k) return(character(0))
    unique(substring(s, seq_len(n - k + 1), seq_len(n - k + 1) + k - 1))
  })
}

#' Test motif enrichment against background via Fisher / binomial test
#' @keywords internal
.test_motif_enrichment <- function(candidates, kmer_map, seq_ids,
                                   bg, n_input, k, chain_label) {

  results <- vector("list", length(candidates))

  pb <- .bluster_progress_new(length(candidates),
                              getOption("bluster.verbose", TRUE),
                              label = sprintf(
                                "[blusteR] Testing %d exact k=%d motifs (%s)...",
                                length(candidates), k, chain_label))
  pb_i <- 0L

  for (i in seq_along(candidates)) {
    motif <- candidates[i]

    # Count input sequences containing this motif
    has_motif <- vapply(kmer_map, function(x) motif %in% x, logical(1))
    n_obs <- sum(has_motif)
    freq_obs <- n_obs / n_input

    # Background frequency
    freq_ref <- bg$per_seq[motif]
    if (is.na(freq_ref)) freq_ref <- 1 / bg$n_total  # pseudocount

    # One-sided binomial test for enrichment
    pval <- stats::binom.test(
      x = n_obs,
      n = n_input,
      p = freq_ref,
      alternative = "greater"
    )$p.value

    fold <- freq_obs / max(freq_ref, 1e-10)

    results[[i]] <- data.table::data.table(
      motif           = motif,
      k               = k,
      chain           = chain_label,
      count_obs       = n_obs,
      freq_obs        = freq_obs,
      freq_ref        = freq_ref,
      fold_enrichment = fold,
      pvalue          = pval,
      pvalue_adj      = NA_real_,  # filled later
      member_ids      = list(seq_ids[has_motif])
    )

    pb_i <- .bluster_progress_tick(pb, pb_i)
  }

  .bluster_progress_close(pb)

  data.table::rbindlist(results)
}

#' Generate degenerate motifs (one position → physico-chemical group)
#' @keywords internal
.generate_degenerate_motifs <- function(motifs, k) {

  # Create group-symbol mapping
  group_symbols <- list()
  for (gname in names(.AA_GROUPS)) {
    sym <- paste0("[", paste0(.AA_GROUPS[[gname]], collapse = ""), "]")
    for (aa in .AA_GROUPS[[gname]]) {
      group_symbols[[aa]] <- sym
    }
  }

  unique_degen <- character(0)

  for (motif in motifs) {
    chars <- strsplit(motif, "")[[1]]
    for (pos in seq_along(chars)) {
      if (chars[pos] %in% names(group_symbols)) {
        new_chars <- chars
        new_chars[pos] <- group_symbols[[chars[pos]]]
        degen <- paste0(new_chars, collapse = "")
        unique_degen <- c(unique_degen, degen)
      }
    }
  }

  unique(unique_degen)
}

#' Extract degenerate motif matches per sequence
#' @keywords internal
.extract_degenerate_kmers_per_seq <- function(cdr3_vec, degen_motifs) {

  pb <- .bluster_progress_new(
    length(cdr3_vec), getOption("bluster.verbose", TRUE),
    label = sprintf("[blusteR] Scanning %d sequences for %d degenerate motifs...",
                    length(cdr3_vec), length(degen_motifs)))
  pb_i <- 0L

  out <- lapply(cdr3_vec, function(s) {
    matches <- vapply(degen_motifs, function(dm) {
      grepl(dm, s)
    }, logical(1))
    pb_i <<- .bluster_progress_tick(pb, pb_i)
    degen_motifs[matches]
  })

  .bluster_progress_close(pb)
  out
}

#' Count input sequences matching each degenerate motif
#' @keywords internal
.count_degenerate_motifs <- function(degen_motifs, degen_map) {

  pb <- .bluster_progress_new(
    length(degen_motifs), getOption("bluster.verbose", TRUE),
    label = sprintf("[blusteR] Counting matches for %d degenerate motifs...",
                    length(degen_motifs)))
  pb_i <- 0L

  counts <- vapply(degen_motifs, function(dm) {
    n <- sum(vapply(degen_map, function(x) dm %in% x, logical(1)))
    pb_i <<- .bluster_progress_tick(pb, pb_i)
    n
  }, integer(1))

  .bluster_progress_close(pb)
  counts
}

#' Test degenerate motif enrichment
#' @keywords internal
.test_degenerate_enrichment <- function(candidates, degen_map, seq_ids,
                                        bg, n_input, k, chain_label) {

  results <- vector("list", length(candidates))

  pb <- .bluster_progress_new(
    length(candidates), getOption("bluster.verbose", TRUE),
    label = sprintf("[blusteR] Testing %d degenerate k=%d motifs (%s)...",
                    length(candidates), k, chain_label))
  pb_i <- 0L

  for (i in seq_along(candidates)) {
    motif <- candidates[i]
    has_motif <- vapply(degen_map, function(x) motif %in% x, logical(1))
    n_obs <- sum(has_motif)
    freq_obs <- n_obs / n_input

    # Real background: fraction of background sequences whose exact k-mers
    # match this degenerate pattern (summed over matching exact k-mers,
    # capped < 1, with a pseudocount floor).
    matching <- grepl(motif, names(bg$per_seq))
    freq_ref <- sum(bg$per_seq[matching])
    if (!is.finite(freq_ref) || freq_ref <= 0) freq_ref <- 1 / bg$n_total
    freq_ref <- min(freq_ref, 0.99)

    pval <- stats::binom.test(
      x = n_obs, n = n_input, p = max(freq_ref, 1e-10),
      alternative = "greater"
    )$p.value

    fold <- freq_obs / max(freq_ref, 1e-10)

    results[[i]] <- data.table::data.table(
      motif           = motif,
      k               = k,
      chain           = chain_label,
      count_obs       = n_obs,
      freq_obs        = freq_obs,
      freq_ref        = freq_ref,
      fold_enrichment = fold,
      pvalue          = pval,
      pvalue_adj      = NA_real_,
      member_ids      = list(seq_ids[has_motif])
    )

    pb_i <- .bluster_progress_tick(pb, pb_i)
  }

  .bluster_progress_close(pb)

  data.table::rbindlist(results)
}
