#' Filter BCR sequences for blusteR analysis
#'
#' Applies quality filters appropriate for specificity-group clustering:
#' removes sequences with ambiguous amino acids, extreme CDR3 lengths,
#' or missing V-gene assignments.
#'
#' @param bcr_data blusteR standard format data.table.
#' @param min_cdr3_len Minimum CDR3 amino acid length (default 5).
#' @param max_cdr3_len Maximum CDR3 amino acid length (default 35).
#' @param require_vgene Require V-gene annotation (default TRUE).
#' @param remove_ambiguous Remove sequences with X or * in CDR3
#'   (default TRUE).
#'
#' @return Filtered data.table.
#' @export
filter_bcr <- function(bcr_data,
                       min_cdr3_len = 5L,
                       max_cdr3_len = 35L,
                       require_vgene = TRUE,
                       remove_ambiguous = TRUE) {

  n_before <- nrow(bcr_data)
  dt <- data.table::copy(bcr_data)

  # CDR3 length filter
  dt <- dt[nchar(cdr3_aa) >= min_cdr3_len & nchar(cdr3_aa) <= max_cdr3_len]

  # Ambiguous residues
  if (remove_ambiguous) {
    dt <- dt[!grepl("[X\\*_]", cdr3_aa)]
  }

  # V-gene requirement
  if (require_vgene) {
    dt <- dt[!is.na(v_gene) & v_gene != "" & v_gene != "None"]
  }

  n_after <- nrow(dt)
  message(sprintf("[blusteR] Filtering: %d → %d contigs (%d removed).",
                  n_before, n_after, n_before - n_after))
  dt[]
}


#' Summarise blusteR clusters
#'
#' Prints a human-readable summary of cluster results.
#'
#' @param bluster_result A \code{bluster_result} object.
#' @param top_n Number of top clusters to show (default 20).
#'
#' @return Invisibly, the clusters data.table.
#' @export
summarize_clusters <- function(bluster_result, top_n = 20L) {

  stopifnot(inherits(bluster_result, "bluster_result"))

  cl <- bluster_result$clusters

  if (nrow(cl) == 0) {
    message("No clusters found.")
    return(invisible(cl))
  }

  cat("\n=== blusteR Cluster Summary ===\n\n")
  cat(sprintf("Total clusters:       %d\n", nrow(cl)))
  cat(sprintf("Total cells assigned: %d\n", nrow(bluster_result$membership)))
  cat(sprintf("Cluster sizes:        %d – %d (median %d)\n",
              min(cl$n_members), max(cl$n_members),
              as.integer(median(cl$n_members))))
  cat(sprintf("Multi-sample:         %d clusters\n",
              sum(cl$n_samples > 1)))

  if (!is.null(bluster_result$annotations) &&
      nrow(bluster_result$annotations) > 0) {
    n_ann <- length(unique(bluster_result$annotations$cluster_id))
    cat(sprintf("Annotated:            %d clusters\n", n_ann))
  }

  cat(sprintf("\nTop %d clusters by blusteR score:\n\n", min(top_n, nrow(cl))))

  show_cols <- intersect(
    c("cluster_id", "n_members", "n_clones", "n_samples",
      "cdr3h_consensus", "top_vh", "bluster_score"),
    names(cl)
  )
  print(head(cl[, ..show_cols], top_n))

  invisible(cl)
}


#' Export blusteR results to files
#'
#' Writes cluster membership, cluster summaries, and annotations to
#' CSV files.
#'
#' @param bluster_result A \code{bluster_result} object.
#' @param prefix Output file prefix (default \code{"bluster_results"}).
#' @param dir Output directory (default \code{"."}).
#'
#' @return Invisibly, a character vector of written file paths.
#' @export
export_clusters <- function(bluster_result,
                            prefix = "bluster_results",
                            dir = ".") {

  stopifnot(inherits(bluster_result, "bluster_result"))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  paths <- character(0)

  # Cluster summaries
  f1 <- file.path(dir, paste0(prefix, "_clusters.csv"))
  cl_export <- data.table::copy(bluster_result$clusters)
  # Remove list columns for CSV
  list_cols <- names(cl_export)[vapply(cl_export, is.list, logical(1))]
  if (length(list_cols) > 0) cl_export[, (list_cols) := NULL]
  utils::write.csv(cl_export, f1, row.names = FALSE)
  paths <- c(paths, f1)

  # Membership
  f2 <- file.path(dir, paste0(prefix, "_membership.csv"))
  utils::write.csv(bluster_result$membership, f2, row.names = FALSE)
  paths <- c(paths, f2)

  # Motifs
  if (nrow(bluster_result$motifs) > 0) {
    f3 <- file.path(dir, paste0(prefix, "_motifs.csv"))
    motifs_export <- data.table::copy(bluster_result$motifs)
    list_cols <- names(motifs_export)[vapply(motifs_export, is.list, logical(1))]
    if (length(list_cols) > 0) {
      for (lc in list_cols) {
        motifs_export[[lc]] <- vapply(motifs_export[[lc]],
                                      function(x) paste(x, collapse = ";"),
                                      character(1))
      }
    }
    utils::write.csv(motifs_export, f3, row.names = FALSE)
    paths <- c(paths, f3)
  }

  # Annotations
  if (!is.null(bluster_result$annotations) &&
      nrow(bluster_result$annotations) > 0) {
    f4 <- file.path(dir, paste0(prefix, "_annotations.csv"))
    utils::write.csv(bluster_result$annotations, f4, row.names = FALSE)
    paths <- c(paths, f4)
  }

  message("[blusteR] Results exported to: ", dir)
  invisible(paths)
}
