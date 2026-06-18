#' Compute global CDR3 similarity edges
#'
#' Analogous to the GLIPH2 global-similarity step, but adapted for BCR
#' biology.  Key differences:
#'
#' \enumerate{
#'   \item \strong{SHM-aware distance}: Instead of raw Hamming distance,
#'     uses BLOSUM62-weighted substitution scoring to distinguish
#'     conservative SHM substitutions (likely shared specificity) from
#'     radical changes (likely different specificity).
#'   \item \strong{Length tolerance}: Allows alignment of CDR3 sequences
#'     differing by up to 1 amino acid in length (BCR CDR3H has more
#'     length variation than TCR CDR3β).
#'   \item \strong{Paired scoring}: When both heavy and light chain CDR3
#'     are available, computes a combined distance.
#' }
#'
#' @param bcr_data A \code{data.table} in blusteR standard format.
#' @param max_dist Maximum SHM-weighted distance for an edge (default 2
#'   for heavy, 1 for light).
#' @param length_diff Maximum allowed CDR3 length difference (default 1).
#' @param chain Which chain: \code{"heavy"} (default), \code{"light"},
#'   or \code{"paired"}.
#' @param scoring \code{"blosum"} (default) for BLOSUM62-weighted
#'   distance, or \code{"hamming"} for unweighted.
#' @param subsample If the number of unique CDR3 sequences exceeds this
#'   value, subsample to reduce computation (default 10000).
#'
#' @return A \code{data.table} of edges with columns: \code{seq_id_1,
#'   seq_id_2, cdr3_1, cdr3_2, distance, chain}.
#'
#' @export
bluster_global <- function(bcr_data,
                         max_dist = .DEFAULT_GLOBAL_DIST,
                         length_diff = 1L,
                         chain = c("heavy", "light", "paired"),
                         scoring = c("blosum", "hamming"),
                         subsample = 10000L) {

  chain <- match.arg(chain)
  scoring <- match.arg(scoring)

  if (chain == "paired") {
    return(.global_paired(bcr_data, max_dist, length_diff, scoring, subsample))
  }

  # Select the appropriate chain
  ch_type <- if (chain == "heavy") "heavy" else "light"
  seqs <- bcr_data[chain_type == ch_type, .(
    seq_id = cell_id,
    cdr3   = cdr3_aa,
    v_gene = v_gene
  )]
  seqs <- unique(seqs, by = "cdr3")
  seqs <- seqs[!is.na(cdr3) & nchar(cdr3) >= 5]

  if (nrow(seqs) == 0) {
    message("[blusteR] No ", chain, " chain sequences for global comparison.")
    return(.empty_edges())
  }

  if (nrow(seqs) > subsample) {
    message(sprintf("[blusteR] Subsampling from %d to %d unique CDR3s.",
                    nrow(seqs), subsample))
    seqs <- seqs[sample(.N, subsample)]
  }

  message(sprintf("[blusteR] Computing global %s similarity for %d unique CDR3s...",
                  chain, nrow(seqs)))

  # Group by CDR3 length (± length_diff)
  seqs[, cdr3_len := nchar(cdr3)]
  edges <- .compute_global_edges(seqs, max_dist, length_diff, scoring, chain)

  message(sprintf("[blusteR] Found %d global similarity edges.", nrow(edges)))
  edges[]
}


# ---- internal global-similarity helpers ------------------------------

#' Compute edges between CDR3 sequences within distance threshold
#' @keywords internal
.compute_global_edges <- function(seqs, max_dist, length_diff, scoring, chain_label) {

  edges_list <- list()
  lengths <- sort(unique(seqs$cdr3_len))

  for (l in lengths) {
    # Get sequences of this length and nearby lengths
    allowed_lens <- seq(l, l + length_diff)
    subset_idx <- which(seqs$cdr3_len %in% allowed_lens)
    if (length(subset_idx) < 2) next

    sub <- seqs[subset_idx]

    if (scoring == "blosum") {
      edge_dt <- .blosum_edges(sub, max_dist, chain_label)
    } else {
      edge_dt <- .hamming_edges(sub, max_dist, chain_label)
    }

    if (nrow(edge_dt) > 0) {
      edges_list[[length(edges_list) + 1]] <- edge_dt
    }
  }

  if (length(edges_list) == 0) return(.empty_edges())
  data.table::rbindlist(edges_list)
}

#' Compute BLOSUM62-weighted distance between CDR3 pairs
#' @keywords internal
.blosum_edges <- function(seqs, max_dist, chain_label) {

  blosum <- .get_blosum62()
  n <- nrow(seqs)
  edges <- list()

  # Pairwise comparison
  for (i in seq_len(n - 1)) {
    s1 <- strsplit(seqs$cdr3[i], "")[[1]]
    l1 <- length(s1)

    for (j in (i + 1):n) {
      s2 <- strsplit(seqs$cdr3[j], "")[[1]]
      l2 <- length(s2)

      # Align sequences (handle length diff by simple end-gap)
      if (l1 == l2) {
        d <- .blosum_distance(s1, s2, blosum)
      } else {
        d <- .blosum_distance_gapped(s1, s2, blosum)
      }

      if (d <= max_dist) {
        edges[[length(edges) + 1]] <- data.table::data.table(
          seq_id_1 = seqs$seq_id[i],
          seq_id_2 = seqs$seq_id[j],
          cdr3_1   = seqs$cdr3[i],
          cdr3_2   = seqs$cdr3[j],
          distance = d,
          chain    = chain_label
        )
      }
    }
  }

  if (length(edges) == 0) return(.empty_edges())
  data.table::rbindlist(edges)
}

#' BLOSUM62 distance: count positions where BLOSUM62 score < 0
#' @keywords internal
.blosum_distance <- function(s1, s2, blosum) {
  min_len <- min(length(s1), length(s2))
  s1 <- s1[seq_len(min_len)]
  s2 <- s2[seq_len(min_len)]

  mismatches <- 0L
  for (pos in seq_len(min_len)) {
    a1 <- s1[pos]
    a2 <- s2[pos]
    if (a1 == a2) next
    # A substitution "costs" 1 if BLOSUM62 score is negative (non-conservative)
    # and 0.5 if score is 0 (neutral), 0 if positive (conservative/SHM-like)
    score <- tryCatch(blosum[a1, a2], error = function(e) -1)
    if (score < 0) {
      mismatches <- mismatches + 1L
    } else if (score == 0) {
      mismatches <- mismatches + 0.5
    }
    # Positive BLOSUM62 score → conservative substitution, don't penalise
  }
  mismatches
}

#' BLOSUM62 distance with one-position gap (indel)
#' @keywords internal
.blosum_distance_gapped <- function(s1, s2, blosum) {

  # Try all possible single-gap positions, take minimum distance
  if (length(s1) > length(s2)) {
    short <- s2; long <- s1
  } else {
    short <- s1; long <- s2
  }

  best <- Inf
  ls <- length(short)
  ll <- length(long)

  if (ll - ls > 1) return(Inf)  # too different

  # Try each gap position in the longer sequence
  for (gap_pos in seq_len(ll)) {
    aligned_long <- long[-gap_pos]
    d <- .blosum_distance(short, aligned_long, blosum) + 0.5  # gap penalty
    best <- min(best, d)
  }

  best
}

#' Hamming distance edges (simpler, faster)
#' @keywords internal
.hamming_edges <- function(seqs, max_dist, chain_label) {

  # Group by exact length for Hamming
  edges <- list()
  for (l in unique(seqs$cdr3_len)) {
    sub <- seqs[cdr3_len == l]
    if (nrow(sub) < 2) next

    mat <- stringdist::stringdistmatrix(sub$cdr3, sub$cdr3, method = "hamming")
    idx <- which(mat <= max_dist & mat > 0, arr.ind = TRUE)
    idx <- idx[idx[,1] < idx[,2], , drop = FALSE]

    if (nrow(idx) > 0) {
      for (r in seq_len(nrow(idx))) {
        edges[[length(edges) + 1]] <- data.table::data.table(
          seq_id_1 = sub$seq_id[idx[r, 1]],
          seq_id_2 = sub$seq_id[idx[r, 2]],
          cdr3_1   = sub$cdr3[idx[r, 1]],
          cdr3_2   = sub$cdr3[idx[r, 2]],
          distance = mat[idx[r, 1], idx[r, 2]],
          chain    = chain_label
        )
      }
    }
  }

  if (length(edges) == 0) return(.empty_edges())
  data.table::rbindlist(edges)
}

#' Paired heavy+light chain global similarity
#' @keywords internal
.global_paired <- function(bcr_data, max_dist, length_diff, scoring, subsample) {

  # Build paired table
  heavy <- bcr_data[chain_type == "heavy", .(
    cell_id, cdr3h = cdr3_aa, vh = v_gene)]
  light <- bcr_data[chain_type == "light", .(
    cell_id, cdr3l = cdr3_aa, vl = v_gene)]

  paired <- merge(heavy, light, by = "cell_id")
  paired <- paired[!is.na(cdr3h) & !is.na(cdr3l)]

  if (nrow(paired) < 2) {
    message("[blusteR] Fewer than 2 paired sequences; cannot do paired global.")
    return(.empty_edges())
  }

  if (nrow(paired) > subsample) {
    paired <- paired[sample(.N, subsample)]
  }

  message(sprintf("[blusteR] Computing paired global similarity for %d cells...",
                  nrow(paired)))

  blosum <- .get_blosum62()
  edges <- list()
  n <- nrow(paired)

  for (i in seq_len(n - 1)) {
    h1 <- strsplit(paired$cdr3h[i], "")[[1]]
    l1 <- strsplit(paired$cdr3l[i], "")[[1]]

    for (j in (i + 1):n) {
      h2 <- strsplit(paired$cdr3h[j], "")[[1]]
      l2 <- strsplit(paired$cdr3l[j], "")[[1]]

      # Combined distance: heavy + light chain
      dh <- if (abs(length(h1) - length(h2)) <= length_diff) {
        if (length(h1) == length(h2)) .blosum_distance(h1, h2, blosum)
        else .blosum_distance_gapped(h1, h2, blosum)
      } else Inf

      dl <- if (abs(length(l1) - length(l2)) <= length_diff) {
        if (length(l1) == length(l2)) .blosum_distance(l1, l2, blosum)
        else .blosum_distance_gapped(l1, l2, blosum)
      } else Inf

      # Combined score: weighted sum (heavy contributes more to specificity)
      d_combined <- 0.7 * dh + 0.3 * dl

      if (d_combined <= max_dist) {
        edges[[length(edges) + 1]] <- data.table::data.table(
          seq_id_1 = paired$cell_id[i],
          seq_id_2 = paired$cell_id[j],
          cdr3_1   = paste0(paired$cdr3h[i], "|", paired$cdr3l[i]),
          cdr3_2   = paste0(paired$cdr3h[j], "|", paired$cdr3l[j]),
          distance = d_combined,
          chain    = "paired"
        )
      }
    }
  }

  if (length(edges) == 0) return(.empty_edges())
  out <- data.table::rbindlist(edges)
  message(sprintf("[blusteR] Found %d paired global edges.", nrow(out)))
  out[]
}

#' Empty edge table constructor
#' @keywords internal
.empty_edges <- function() {
  data.table::data.table(
    seq_id_1 = character(0), seq_id_2 = character(0),
    cdr3_1 = character(0), cdr3_2 = character(0),
    distance = numeric(0), chain = character(0)
  )
}
