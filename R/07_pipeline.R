#' Run the complete blusteR pipeline
#'
#' One-call convenience function that runs local motif discovery,
#' global similarity, clustering, scoring, and optional epitope
#' annotation.
#'
#' @param bcr_data A \code{data.table} in blusteR standard format (output
#'   of \code{load_10x_bcr()}, \code{load_screpertoire()}, or
#'   \code{load_airr()}).
#' @param reference A \code{bluster_reference} object.  If \code{NULL},
#'   the pre-built reference bundled with the package is loaded via
#'   \code{load_reference()}.
#' @param chain Which chain to analyse: \code{"heavy"} (default),
#'   \code{"light"}, \code{"both"}, or \code{"paired"}.
#' @param k_sizes K-mer lengths for motif discovery (default
#'   \code{c(4, 5)}).
#' @param motif_min_freq Minimum motif frequency (default 3).
#' @param motif_max_freq Maximum fraction of input sequences a motif may
#'   appear in (default 0.20); more common motifs are dropped.
#' @param motif_p_cutoff Adjusted p-value threshold for motifs
#'   (default 0.05).
#' @param global_max_dist Maximum global CDR3 distance (default 2).
#' @param global_scoring \code{"blosum"} or \code{"hamming"}.
#' @param min_cluster_size Minimum cluster members (default 3).
#' @param clustering_method \code{"components"} or \code{"louvain"}.
#' @param collapse_clones Collapse clonal relatives first (default TRUE).
#' @param annotate Logical; annotate against IEDB/SAbDab (default TRUE).
#' @param species \code{"human"} or \code{"mouse"}.
#' @param verbose Print progress messages (default TRUE).
#'
#' @return A \code{bluster_result} object.
#'
#' @examples
#' \dontrun{
#' # From 10X Genomics VDJ output
#' bcr <- load_10x_bcr("/path/to/cellranger/vdj/output/")
#' result <- bluster(bcr)
#' print(result$clusters)
#' plot_cluster_network(result)
#'
#' # From scRepertoire
#' library(scRepertoire)
#' combined <- combineBCR(contig_list, samples = c("S1", "S2"))
#' bcr <- load_screpertoire(combined)
#' result <- bluster(bcr, chain = "paired")
#'
#' # From multiple 10X samples
#' bcr1 <- load_10x_bcr("sample1/", sample_id = "donor1")
#' bcr2 <- load_10x_bcr("sample2/", sample_id = "donor2")
#' bcr <- rbind(bcr1, bcr2)
#' result <- bluster(bcr, chain = "both", annotate = TRUE)
#' }
#'
#' @export
bluster <- function(bcr_data,
                  reference = NULL,
                  chain = c("heavy", "light", "both", "paired"),
                  k_sizes = .DEFAULT_KMER_SIZES,
                  motif_min_freq = .DEFAULT_LOCAL_MIN_FREQ,
                  motif_max_freq = .DEFAULT_LOCAL_MAX_FREQ,
                  motif_p_cutoff = .DEFAULT_PVALUE,
                  global_max_dist = .DEFAULT_GLOBAL_DIST,
                  global_scoring = c("blosum", "hamming"),
                  min_cluster_size = 3L,
                  clustering_method = c("components", "louvain"),
                  collapse_clones = TRUE,
                  annotate = TRUE,
                  species = c("human", "mouse"),
                  verbose = TRUE) {

  chain <- match.arg(chain)
  global_scoring <- match.arg(global_scoring)
  clustering_method <- match.arg(clustering_method)
  species <- match.arg(species)

  if (!verbose) {
    old_opts <- options(bluster.verbose = FALSE)
    on.exit(options(old_opts), add = TRUE)
  }

  # Validate input
  .validate_bcr_data(bcr_data)

  message("╔══════════════════════════════════════════════════════════╗")
  message("║            blusteR — B-cell Interaction by               ║")
  message("║                Paratope Hotspots v1.0                    ║")
  message("╚══════════════════════════════════════════════════════════╝")
  message("")

  # --- Step 0: Build reference ------------------------------------------
  if (is.null(reference)) {
    message("── Step 0: Loading reference databases ──────────────────")
    reference <- load_reference(species = species)
    message("")
  }

  # --- Step 1: Pre-filter -----------------------------------------------
  message("── Step 1: Pre-filtering BCR data ────────────────────────")
  bcr_filtered <- filter_bcr(bcr_data)
  message("")

  # --- Step 2: Local motif discovery ------------------------------------
  message("── Step 2: Local motif discovery ─────────────────────────")
  motif_chain <- if (chain == "paired") "both" else chain
  motifs <- bluster_motifs(
    bcr_data     = bcr_filtered,
    reference    = reference,
    chain        = motif_chain,
    k_sizes      = k_sizes,
    min_freq     = motif_min_freq,
    max_freq     = motif_max_freq,
    p_cutoff     = motif_p_cutoff
  )
  message("")

  # --- Step 3: Global similarity ----------------------------------------
  message("── Step 3: Global CDR3 similarity ────────────────────────")
  global_chain <- if (chain == "both") "heavy" else chain
  global <- bluster_global(
    bcr_data  = bcr_filtered,
    max_dist  = global_max_dist,
    chain     = global_chain,
    scoring   = global_scoring
  )
  message("")

  # --- Step 4: Clustering -----------------------------------------------
  message("── Step 4: Cluster formation & scoring ───────────────────")
  result <- bluster_cluster(
    bcr_data          = bcr_filtered,
    motif_results     = motifs,
    global_edges      = global,
    min_cluster_size  = min_cluster_size,
    clustering_method = clustering_method,
    collapse_clones   = collapse_clones
  )
  message("")

  # --- Step 5: Epitope annotation (optional) ----------------------------
  if (annotate && nrow(result$clusters) > 0) {
    message("── Step 5: Epitope annotation ───────────────────────────")
    result <- annotate_epitopes(
      bluster_result = result,
      reference    = reference,
      bcr_data     = bcr_filtered
    )
    message("")
  }

  # --- Summary ----------------------------------------------------------
  message("══════════════════════════════════════════════════════════")
  message(sprintf("  Input:      %d cells, %d samples",
                  length(unique(bcr_filtered$cell_id)),
                  length(unique(bcr_filtered$sample_id))))
  message(sprintf("  Motifs:     %d enriched", nrow(motifs)))
  message(sprintf("  Global:     %d edges", nrow(global)))
  message(sprintf("  Clusters:   %d (≥%d members)",
                  nrow(result$clusters), min_cluster_size))
  if (!is.null(result$annotations) && nrow(result$annotations) > 0) {
    message(sprintf("  Annotated:  %d clusters with known specificity",
                    length(unique(result$annotations$cluster_id))))
  }
  message("══════════════════════════════════════════════════════════")

  result
}


#' Create a pipeline progress bar
#'
#' Returns a base-R text progress bar tracking progress over \code{total}
#' iterations, or \code{NULL} when \code{verbose} is \code{FALSE}.  An
#' optional \code{label} is printed on its own line before the bar so the
#' user can tell which section is running.
#' @keywords internal
.bluster_progress_new <- function(total, verbose = TRUE, label = NULL) {
  if (!isTRUE(verbose) || total < 1L) return(NULL)
  if (!is.null(label)) message(label)
  utils::txtProgressBar(min = 0, max = total, initial = 0, style = 3)
}

#' Advance a pipeline progress bar by one step
#'
#' @return The new (incremented) step count.
#' @keywords internal
.bluster_progress_tick <- function(pb, step_i) {
  step_i <- step_i + 1L
  if (!is.null(pb)) utils::setTxtProgressBar(pb, step_i)
  step_i
}

#' Close a pipeline progress bar
#' @keywords internal
.bluster_progress_close <- function(pb) {
  if (!is.null(pb)) close(pb)
  invisible(NULL)
}


#' Validate blusteR standard-format BCR data
#' @keywords internal
.validate_bcr_data <- function(dt) {

  if (!data.table::is.data.table(dt)) {
    stop("Input must be a data.table. Use load_10x_bcr(), load_screpertoire(), or load_airr().")
  }

  required <- c("cell_id", "chain_type", "cdr3_aa")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "),
         "\nUse load_10x_bcr(), load_screpertoire(), or load_airr() to prepare input.")
  }

  if (nrow(dt) == 0) stop("Input data has 0 rows.")

  n_heavy <- sum(dt$chain_type == "heavy", na.rm = TRUE)
  n_light <- sum(dt$chain_type == "light", na.rm = TRUE)
  message(sprintf("[blusteR] Input: %d heavy chain, %d light chain contigs.",
                  n_heavy, n_light))

  invisible(TRUE)
}
