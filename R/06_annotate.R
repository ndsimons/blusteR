#' Annotate blusteR clusters with known epitope specificities
#'
#' Cross-references cluster CDR3 sequences against IEDB B-cell epitope
#' records and SAbDab structural data to identify clusters whose members
#' resemble antibodies of known specificity.
#'
#' Matching is performed by:
#' \enumerate{
#'   \item \strong{Exact CDR3 match} against IEDB/SAbDab CDR3H sequences.
#'   \item \strong{Near match} (BLOSUM62 distance ≤ threshold) for CDR3H.
#'   \item \strong{Paired match} using both CDR3H and CDR3L when available.
#' }
#'
#' @param bluster_result Output from \code{bluster_cluster()}.
#' @param reference A \code{bluster_reference} object.
#' @param max_dist Maximum BLOSUM62 distance for near-match annotation
#'   (default 1).
#' @param bcr_data Original BCR data for light-chain lookup.
#'
#' @return The input \code{bluster_result} with an added
#'   \code{$annotations} data.table containing cluster-level antigen
#'   annotations.
#'
#' @export
annotate_epitopes <- function(bluster_result,
                              reference = NULL,
                              max_dist = 1L,
                              bcr_data = NULL) {

  stopifnot(inherits(bluster_result, "bluster_result"))

  if (is.null(reference)) {
    message("[blusteR] Building reference for annotation...")
    reference <- build_reference()
  }

  clusters <- bluster_result$clusters
  membership <- bluster_result$membership

  if (nrow(clusters) == 0) {
    bluster_result$annotations <- data.table::data.table()
    return(bluster_result)
  }

  message("[blusteR] Annotating clusters against IEDB and SAbDab...")

  # Collect CDR3H per cluster
  if (!is.null(bcr_data)) {
    cluster_data <- merge(membership, bcr_data, by = "cell_id", all.x = TRUE)
  } else {
    cluster_data <- membership
  }

  annotations <- list()

  # Build combined known-specificity table from IEDB + SAbDab
  known <- .build_known_specificity_table(reference)

  if (nrow(known) == 0) {
    message("[blusteR] No known specificities in reference to match against.")
    bluster_result$annotations <- data.table::data.table()
    return(bluster_result)
  }

  blosum <- .get_blosum62()

  for (cl_id in unique(membership$cluster_id)) {
    cl_cells <- membership[cluster_id == cl_id, cell_id]

    if (!is.null(bcr_data)) {
      cl_cdr3h <- unique(bcr_data[cell_id %in% cl_cells &
                                    chain_type == "heavy", cdr3_aa])
    } else {
      cl_cdr3h <- clusters[cluster_id == cl_id, cdr3h_consensus]
    }
    cl_cdr3h <- cl_cdr3h[!is.na(cl_cdr3h)]

    for (cdr3 in cl_cdr3h) {
      cdr3_canon <- .canon_cdr3(cdr3)
      if (is.na(cdr3_canon) || nchar(cdr3_canon) == 0) next
      s1 <- strsplit(cdr3_canon, "")[[1]]

      # Exact match (boundary-normalised)
      exact <- known[cdr3h_canon == cdr3_canon]
      if (nrow(exact) > 0) {
        for (r in seq_len(nrow(exact))) {
          annotations[[length(annotations) + 1]] <- data.table::data.table(
            cluster_id = cl_id,
            query_cdr3 = cdr3,
            match_cdr3 = exact$cdr3h[r],
            match_type = "exact",
            distance   = 0,
            antigen    = exact$antigen[r],
            organism   = exact$organism[r],
            source     = exact$source[r],
            pdb_id     = exact$pdb_id[r]
          )
        }
        next
      }

      # Near match (boundary-normalised)
      for (r in seq_len(nrow(known))) {
        ref_canon <- known$cdr3h_canon[r]
        if (is.na(ref_canon) || nchar(ref_canon) == 0) next
        if (abs(nchar(cdr3_canon) - nchar(ref_canon)) > 1) next

        s2 <- strsplit(ref_canon, "")[[1]]
        d <- if (length(s1) == length(s2)) {
          .blosum_distance(s1, s2, blosum)
        } else {
          .blosum_distance_gapped(s1, s2, blosum)
        }

        if (d <= max_dist) {
          annotations[[length(annotations) + 1]] <- data.table::data.table(
            cluster_id = cl_id,
            query_cdr3 = cdr3,
            match_cdr3 = known$cdr3h[r],
            match_type = "near",
            distance   = d,
            antigen    = known$antigen[r],
            organism   = known$organism[r],
            source     = known$source[r],
            pdb_id     = known$pdb_id[r]
          )
        }
      }
    }
  }

  if (length(annotations) > 0) {
    ann_dt <- data.table::rbindlist(annotations, fill = TRUE)
    message(sprintf("[blusteR] Annotated %d clusters with %d known specificities.",
                    length(unique(ann_dt$cluster_id)), nrow(ann_dt)))
  } else {
    ann_dt <- data.table::data.table()
    message("[blusteR] No clusters matched known specificities.")
  }

  bluster_result$annotations <- ann_dt
  bluster_result
}


#' Score existing clusters (convenience wrapper)
#'
#' Re-scores clusters, e.g. after filtering. Runs V-gene enrichment
#' and recomputes composite scores.
#'
#' @param bluster_result A \code{bluster_result} object.
#' @param bcr_data BCR data for V-gene context.
#' @return Updated \code{bluster_result}.
#' @export
score_clusters <- function(bluster_result, bcr_data) {

  stopifnot(inherits(bluster_result, "bluster_result"))
  # Recompute scores using the same logic as bluster_cluster
  # (This is a convenience re-export; the scores are computed in bluster_cluster)
  bluster_result$clusters[, bluster_score := .compute_bluster_score(
    n_members, clonal_diversity, sample_diversity,
    vh_enrichment_p, cdr3h_len_sd
  )]
  data.table::setorder(bluster_result$clusters, -bluster_score)
  bluster_result
}


# ---- internal annotation helpers ------------------------------------

#' Combine IEDB and SAbDab into a unified known-specificity table
#' @keywords internal
.build_known_specificity_table <- function(reference) {

  rows <- list()

  # IEDB entries with CDR3H
  if (!is.null(reference$iedb) && nrow(reference$iedb) > 0) {
    iedb <- reference$iedb
    if ("cdr3h" %in% names(iedb)) {
      iedb_valid <- iedb[!is.na(cdr3h) & nchar(cdr3h) > 0]
      if (nrow(iedb_valid) > 0) {
        rows[[1]] <- data.table::data.table(
          cdr3h    = iedb_valid$cdr3h,
          cdr3l    = if ("cdr3l" %in% names(iedb_valid)) iedb_valid$cdr3l else NA,
          antigen  = iedb_valid$antigen,
          organism = if ("organism" %in% names(iedb_valid)) iedb_valid$organism else NA,
          source   = "IEDB",
          pdb_id   = if ("pdb_id" %in% names(iedb_valid)) iedb_valid$pdb_id else NA
        )
      }
    }
  }

  # SAbDab entries
  if (!is.null(reference$sabdab) && nrow(reference$sabdab) > 0) {
    sab <- reference$sabdab
    if ("cdr3h" %in% names(sab)) {
      sab_valid <- sab[!is.na(cdr3h) & nchar(cdr3h) > 0]
      if (nrow(sab_valid) > 0) {
        rows[[length(rows) + 1]] <- data.table::data.table(
          cdr3h    = sab_valid$cdr3h,
          cdr3l    = if ("cdr3l" %in% names(sab_valid)) sab_valid$cdr3l else NA,
          antigen  = sab_valid$antigen,
          organism = NA_character_,
          source   = "SAbDab",
          pdb_id   = sab_valid$pdb_id
        )
      }
    }
  }

  if (length(rows) == 0) {
    return(data.table::data.table(
      cdr3h = character(0), cdr3l = character(0),
      antigen = character(0), organism = character(0),
      source = character(0), pdb_id = character(0),
      cdr3h_canon = character(0)
    ))
  }

  tab <- data.table::rbindlist(rows, fill = TRUE)

  # Some sources store several CDR3H joined by ";"; expand to one per row
  tab <- tab[!is.na(cdr3h) & nchar(cdr3h) > 0]
  if (nrow(tab) > 0 && any(grepl(";", tab$cdr3h, fixed = TRUE))) {
    tab <- tab[, .(cdr3h = trimws(unlist(strsplit(cdr3h, ";", fixed = TRUE)))),
               by = setdiff(names(tab), "cdr3h")]
    tab <- tab[nchar(cdr3h) > 0]
  }

  # Boundary-normalised CDR3H so input of any convention (with/without the
  # conserved flanking C...W) matches the reference loop residues.
  tab[, cdr3h_canon := .canon_cdr3(cdr3h)]
  tab[!is.na(cdr3h_canon) & nchar(cdr3h_canon) > 0]
}

#' Canonicalise a CDR3 amino-acid sequence for cross-source matching
#'
#' Databases differ in whether they include the conserved residues that
#' flank the CDR3 loop (a Cys before, a Trp or Phe after).  Stripping a
#' single leading \code{C} and trailing \code{W}/\code{F} lets input
#' (e.g. 10X \code{"CARDRW"}) and IEDB (\code{"ARDR"}) conventions be
#' compared on the same loop residues.
#' @keywords internal
.canon_cdr3 <- function(x) {
  x <- toupper(trimws(x))
  x <- sub("^C", "", x)
  x <- sub("[WF]$", "", x)
  x
}
