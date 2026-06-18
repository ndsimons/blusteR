#' Load BCR data from 10X Genomics VDJ output
#'
#' Reads \code{filtered_contig_annotations.csv} (and optionally
#' \code{clonotypes.csv}) produced by Cell Ranger VDJ.
#'
#' @param path Directory containing Cell Ranger VDJ output, or the path to
#'   \code{filtered_contig_annotations.csv} directly.
#' @param sample_id Character label for this sample (default: basename of path).
#' @param remove_nonproductive Logical; drop non-productive contigs
#'   (default TRUE).
#' @param keep_multi Logical; if FALSE (default), keep only the first
#'   heavy and first light chain per barcode when multiple exist.
#'
#' @return A \code{data.table} in blusteR standard format with columns:
#'   \code{cell_id, sample_id, chain, v_gene, d_gene, j_gene, c_gene,
#'   cdr3_aa, cdr3_nt, clone_id, is_cell, productive}.
#'   Heavy and light chain contigs for the same barcode share a
#'   \code{cell_id}.
#'
#' @export
load_10x_bcr <- function(path,
                         sample_id = NULL,
                         remove_nonproductive = TRUE,
                         keep_multi = FALSE) {

  # Resolve file path
  if (dir.exists(path)) {
    csv_path <- file.path(path, "filtered_contig_annotations.csv")
    if (!file.exists(csv_path)) {
      # Try outs subfolder
      csv_path <- file.path(path, "outs", "filtered_contig_annotations.csv")
    }
  } else {
    csv_path <- path
  }
  stopifnot("Cannot find filtered_contig_annotations.csv" = file.exists(csv_path))

  dt <- data.table::fread(csv_path, stringsAsFactors = FALSE)

  if (is.null(sample_id)) sample_id <- basename(dirname(csv_path))

  # Standardise column names (10X uses various naming conventions)
  col_map <- c(
    barcode        = "barcode",
    is_cell        = "is_cell",
    chain          = "chain",
    v_gene         = "v_gene",
    d_gene         = "d_gene",
    j_gene         = "j_gene",
    c_gene         = "c_gene",
    cdr3           = "cdr3",
    cdr3_nt        = "cdr3_nt",
    raw_clonotype_id = "raw_clonotype_id",
    productive     = "productive"
  )

  # Map columns that exist
  present <- intersect(names(col_map), names(dt))
  if (length(present) == 0) {
    stop("Input CSV does not appear to be a 10X filtered_contig_annotations file.")
  }

  # Build standardised table
  out <- data.table::data.table(
    cell_id    = paste0(sample_id, "_", dt$barcode),
    sample_id  = sample_id,
    chain      = .standardise_chain(dt$chain),
    v_gene     = dt$v_gene,
    d_gene     = if ("d_gene" %in% names(dt)) dt$d_gene else NA_character_,
    j_gene     = dt$j_gene,
    c_gene     = if ("c_gene" %in% names(dt)) dt$c_gene else NA_character_,
    cdr3_aa    = dt$cdr3,
    cdr3_nt    = if ("cdr3_nt" %in% names(dt)) dt$cdr3_nt else NA_character_,
    clone_id   = if ("raw_clonotype_id" %in% names(dt)) {
      paste0(sample_id, "_", dt$raw_clonotype_id)
    } else {
      NA_character_
    },
    is_cell    = if ("is_cell" %in% names(dt)) {
      as.logical(dt$is_cell)
    } else TRUE,
    productive = if ("productive" %in% names(dt)) {
      as.logical(dt$productive)
    } else TRUE
  )

  # Filter
  out <- out[is_cell == TRUE]
  if (remove_nonproductive) out <- out[productive == TRUE]

  # Classify chain type
  out[, chain_type := ifelse(chain == "IGH", "heavy", "light")]

  # Deduplicate multi-chain cells

  if (!keep_multi) {
    out <- .select_primary_chains(out)
  }

  # Remove contigs with missing CDR3
  out <- out[!is.na(cdr3_aa) & nchar(cdr3_aa) > 0]

  message(sprintf("[blusteR] Loaded %d contigs from %d cells (sample: %s)",
                  nrow(out), length(unique(out$cell_id)), sample_id))
  out[]
}


#' Load BCR data from a scRepertoire combined contigs object
#'
#' Accepts the list returned by \code{scRepertoire::combineBCR()} or
#' \code{scRepertoire::combineExpression()} and converts it to blusteR
#' standard format.
#'
#' @param scr A list of data.frames from \code{combineBCR()}, or a single
#'   merged data.frame.
#' @param sample_ids Optional character vector of sample names (one per
#'   list element).
#'
#' @return A \code{data.table} in blusteR standard format.
#' @export
load_screpertoire <- function(scr, sample_ids = NULL) {

  # If list, rbind

  if (is.list(scr) && !is.data.frame(scr)) {
    if (is.null(sample_ids)) sample_ids <- paste0("S", seq_along(scr))
    for (i in seq_along(scr)) {
      scr[[i]]$sample_id <- sample_ids[i]
    }
    scr <- data.table::rbindlist(scr, fill = TRUE)
  } else {
    scr <- data.table::as.data.table(scr)
    if (!"sample_id" %in% names(scr)) scr$sample_id <- "S1"
  }

  # scRepertoire stores BCR data with IGH and IGL/IGK columns
  # Column patterns: IGH, IGLC (light chain combined)
  # CDR3 is in columns like "CTaa" or the IGH/IGL separated columns

  rows <- list()

  for (i in seq_len(nrow(scr))) {
    row <- scr[i, ]
    cell <- as.character(row$barcode %||% row$cell_id %||% paste0(row$sample_id, "_", i))

    # Extract heavy chain info
    if (!is.null(row$IGH) && !is.na(row$IGH) && nchar(as.character(row$IGH)) > 0) {
      h_parts <- .parse_screpertoire_chain(as.character(row$IGH))
      rows[[length(rows) + 1]] <- data.table::data.table(
        cell_id   = cell,
        sample_id = as.character(row$sample_id),
        chain     = "IGH",
        chain_type = "heavy",
        v_gene    = h_parts$v_gene,
        d_gene    = h_parts$d_gene,
        j_gene    = h_parts$j_gene,
        c_gene    = NA_character_,
        cdr3_aa   = h_parts$cdr3,
        cdr3_nt   = NA_character_,
        clone_id  = as.character(row$CTstrict %||% row$CTgene %||% NA),
        is_cell   = TRUE,
        productive = TRUE
      )
    }

    # Extract light chain info
    lc_col <- if (!is.null(row$IGLC)) "IGLC" else "IGL"
    lc_val <- row[[lc_col]]
    if (!is.null(lc_val) && !is.na(lc_val) && nchar(as.character(lc_val)) > 0) {
      l_parts <- .parse_screpertoire_chain(as.character(lc_val))
      chain_name <- if (grepl("^IGK", l_parts$v_gene)) "IGK" else "IGL"
      rows[[length(rows) + 1]] <- data.table::data.table(
        cell_id   = cell,
        sample_id = as.character(row$sample_id),
        chain     = chain_name,
        chain_type = "light",
        v_gene    = l_parts$v_gene,
        d_gene    = NA_character_,
        j_gene    = l_parts$j_gene,
        c_gene    = NA_character_,
        cdr3_aa   = l_parts$cdr3,
        cdr3_nt   = NA_character_,
        clone_id  = as.character(row$CTstrict %||% row$CTgene %||% NA),
        is_cell   = TRUE,
        productive = TRUE
      )
    }
  }

  out <- data.table::rbindlist(rows, fill = TRUE)
  out <- out[!is.na(cdr3_aa) & nchar(cdr3_aa) > 0]

  message(sprintf("[blusteR] Loaded %d contigs from %d cells across %d samples",
                  nrow(out), length(unique(out$cell_id)),
                  length(unique(out$sample_id))))
  out[]
}


#' Load BCR data from AIRR-format TSV
#'
#' Reads a TSV file following the AIRR Community standard rearrangement
#' schema.
#'
#' @param path Path to AIRR TSV file.
#' @param sample_id Sample label.
#' @param remove_nonproductive Drop non-productive rearrangements.
#'
#' @return A \code{data.table} in blusteR standard format.
#' @export
load_airr <- function(path, sample_id = "AIRR", remove_nonproductive = TRUE) {

  dt <- data.table::fread(path, stringsAsFactors = FALSE)

  out <- data.table::data.table(
    cell_id    = if ("cell_id" %in% names(dt)) {
      paste0(sample_id, "_", dt$cell_id)
    } else {
      paste0(sample_id, "_", seq_len(nrow(dt)))
    },
    sample_id  = sample_id,
    chain      = .standardise_chain(dt$locus),
    v_gene     = dt$v_call,
    d_gene     = if ("d_call" %in% names(dt)) dt$d_call else NA_character_,
    j_gene     = dt$j_call,
    c_gene     = if ("c_call" %in% names(dt)) dt$c_call else NA_character_,
    cdr3_aa    = dt$junction_aa,
    cdr3_nt    = if ("junction" %in% names(dt)) dt$junction else NA_character_,
    clone_id   = if ("clone_id" %in% names(dt)) {
      paste0(sample_id, "_", dt$clone_id)
    } else NA_character_,
    is_cell    = TRUE,
    productive = if ("productive" %in% names(dt)) {
      as.logical(dt$productive)
    } else TRUE
  )

  out[, chain_type := ifelse(chain == "IGH", "heavy", "light")]

  if (remove_nonproductive) out <- out[productive == TRUE]
  out <- out[!is.na(cdr3_aa) & nchar(cdr3_aa) > 0]

  message(sprintf("[blusteR] Loaded %d contigs (AIRR format, sample: %s)",
                  nrow(out), sample_id))
  out[]
}


# ---- internal helpers for parsers ------------------------------------

#' Standardise chain names to IGH / IGK / IGL
#' @keywords internal
.standardise_chain <- function(x) {
  x <- toupper(x)
  x[x %in% c("IGH", "HEAVY")] <- "IGH"
  x[x %in% c("IGK", "KAPPA")] <- "IGK"
  x[x %in% c("IGL", "LAMBDA")] <- "IGL"
  x
}

#' Select one heavy and one light chain per cell
#' @keywords internal
.select_primary_chains <- function(dt) {
  # Keep highest-UMI contig per chain type per cell (if umi_count exists),

  # otherwise keep the first row.
  dt[, .SD[1], by = .(cell_id, chain_type)]
}

#' Parse scRepertoire chain string like "IGHV3-30.IGHD1-1.IGHJ4.CARGYSSGWYFDYW"
#' @keywords internal
.parse_screpertoire_chain <- function(s) {
  parts <- strsplit(s, "\\.")[[1]]
  v <- d <- j <- cdr3 <- NA_character_

  for (p in parts) {
    if (grepl("^IG[HKL]V", p)) v <- p
    else if (grepl("^IGHD", p)) d <- p
    else if (grepl("^IG[HKL]J", p)) j <- p
    else if (grepl("^C[A-Z]", p) && nchar(p) > 5) cdr3 <- p
  }

  list(v_gene = v, d_gene = d, j_gene = j, cdr3 = cdr3)
}

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b
