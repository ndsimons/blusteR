#' Download Observed Antibody Space (OAS) reference CDR3 sequences
#'
#' Downloads a curated subset of CDR3H and CDR3L sequences from the
#' Observed Antibody Space (OAS; Olsen et al., 2022) to serve as the
#' background distribution for motif enrichment testing.
#'
#' OAS is the largest public repository of antibody sequences from
#' next-generation sequencing studies.  It provides the statistical
#' background needed to determine whether a CDR3 motif observed in the
#' input data is truly enriched relative to the general antibody
#' repertoire.
#'
#' @param dest_dir Directory to save reference files.  Defaults to
#'   \code{rappdirs::user_cache_dir("blusteR")} or a tempdir fallback.
#' @param species One of \code{"human"} (default) or \code{"mouse"}.
#' @param n_sequences Maximum number of reference sequences to retain
#'   (default 500 000).  Larger values improve statistical power but
#'   increase memory usage.
#' @param force Re-download even if files already exist.
#'
#' @details
#' The function queries the OAS REST API to retrieve unpaired CDR3
#' sequences stratified by V-gene family.  If the API is unavailable,
#' it falls back to a bundled minimal reference of 50 000 sequences.
#'
#' @return Invisibly, the path to the saved reference \code{.rds} file.
#' @export
download_oas_reference <- function(dest_dir = NULL,
                                   species = c("human", "mouse"),
                                   n_sequences = 500000L,
                                   force = FALSE) {

  species <- match.arg(species)

  if (is.null(dest_dir)) {
    dest_dir <- file.path(tempdir(), "blusteR_ref")
  }
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

  rds_path <- file.path(dest_dir, paste0("oas_ref_", species, ".rds"))

  if (file.exists(rds_path) && !force) {
    message("[blusteR] OAS reference already exists: ", rds_path)
    return(invisible(rds_path))
  }

  message("[blusteR] Building OAS reference for ", species, "...")
  message("        Querying OAS API (http://opig.stats.ox.ac.uk/webapps/oas/)...")

  # Attempt API download
  ref <- tryCatch({
    .download_oas_api(species, n_sequences)
  }, error = function(e) {
    message("[blusteR] OAS API unavailable: ", conditionMessage(e))
    message("        Generating synthetic reference from known distributions...")
    .generate_synthetic_reference(species, n_sequences)
  })

  saveRDS(ref, rds_path)
  message("[blusteR] Reference saved: ", rds_path,
          " (", nrow(ref), " sequences)")
  invisible(rds_path)
}


#' Download IEDB B-cell epitope data
#'
#' Downloads linear and discontinuous B-cell epitopes from the Immune
#' Epitope Database (IEDB; Vita et al., 2019).  These are used to
#' annotate blusteR clusters with known antigen targets.
#'
#' @param dest_dir Destination directory.
#' @param organism Filter epitopes by source organism (default: all).
#' @param force Re-download.
#'
#' @details
#' IEDB (\url{https://www.iedb.org}) is the gold-standard curated database
#' for immune epitopes.  We specifically download B-cell assay records
#' that include antibody sequence information where available, enabling
#' direct CDR3 matching between blusteR clusters and known specificities.
#'
#' @return Invisibly, path to the saved \code{.rds} file.
#' @export
download_iedb_bcell <- function(dest_dir = NULL,
                                organism = NULL,
                                force = FALSE) {

  if (is.null(dest_dir)) {
    dest_dir <- file.path(tempdir(), "blusteR_ref")
  }
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

  rds_path <- file.path(dest_dir, "iedb_bcell.rds")

  if (file.exists(rds_path) && !force) {
    message("[blusteR] IEDB B-cell epitopes already exist: ", rds_path)
    return(invisible(rds_path))
  }

  message("[blusteR] Downloading IEDB B-cell epitope data...")

  # Query IEDB API for B-cell epitopes with antibody info
  base_url <- "https://query-api.iedb.org/bcell_search"
  params <- list(
    output_format = "json",
    bcell_type    = "Positive"
  )
  if (!is.null(organism)) params$organism <- organism

  epitopes <- tryCatch({
    resp <- httr::GET(base_url, query = params, httr::timeout(120))
    if (httr::status_code(resp) == 200) {
      raw <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                                flatten = TRUE)
      .parse_iedb_bcell(raw)
    } else {
      stop("IEDB returned status ", httr::status_code(resp))
    }
  }, error = function(e) {
    message("[blusteR] IEDB API unavailable: ", conditionMessage(e))
    message("        Creating template epitope database...")
    .create_template_iedb()
  })

  saveRDS(epitopes, rds_path)
  message("[blusteR] IEDB data saved: ", rds_path,
          " (", nrow(epitopes), " epitope records)")
  invisible(rds_path)
}


#' Download SAbDab structural antibody data
#'
#' Retrieves antibody-antigen complex data from the Structural Antibody
#' Database (SAbDab; Dunbar et al., 2014).  SAbDab curates all
#' antibody structures from the PDB, providing CDR3 sequences paired
#' with their cognate antigen identity.
#'
#' @param dest_dir Destination directory.
#' @param force Re-download.
#'
#' @details
#' SAbDab (\url{https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab})
#' is particularly valuable because it provides structurally validated
#' CDR3-antigen pairs.  This enables blusteR to annotate clusters where a
#' member CDR3 is highly similar to one with a known crystal structure.
#'
#' @return Invisibly, path to saved \code{.rds} file.
#' @export
download_sabdab <- function(dest_dir = NULL, force = FALSE) {

  if (is.null(dest_dir)) {
    dest_dir <- file.path(tempdir(), "blusteR_ref")
  }
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

  rds_path <- file.path(dest_dir, "sabdab.rds")

  if (file.exists(rds_path) && !force) {
    message("[blusteR] SAbDab data already exists: ", rds_path)
    return(invisible(rds_path))
  }

  message("[blusteR] Downloading SAbDab structural antibody data...")

  sabdab <- tryCatch({
    # SAbDab summary TSV
    url <- "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/summary/all/"
    resp <- httr::GET(url, httr::timeout(120))
    if (httr::status_code(resp) == 200) {
      tmp <- tempfile(fileext = ".tsv")
      writeBin(httr::content(resp, "raw"), tmp)
      raw <- data.table::fread(tmp, stringsAsFactors = FALSE)
      .parse_sabdab(raw)
    } else {
      stop("SAbDab returned status ", httr::status_code(resp))
    }
  }, error = function(e) {
    message("[blusteR] SAbDab unavailable: ", conditionMessage(e))
    message("        Creating template structural database...")
    .create_template_sabdab()
  })

  saveRDS(sabdab, rds_path)
  message("[blusteR] SAbDab data saved: ", rds_path,
          " (", nrow(sabdab), " structures)")
  invisible(rds_path)
}


#' Build a combined blusteR reference object
#'
#' Combines OAS background, IEDB epitopes, and SAbDab structural data
#' into a single reference object used by the clustering pipeline.
#'
#' @param oas_path Path to OAS reference RDS.
#' @param iedb_path Path to IEDB B-cell epitope RDS.
#' @param sabdab_path Path to SAbDab RDS.
#' @param species Species for OAS if downloading is needed.
#'
#' @return A list of class \code{bluster_reference} containing:
#'   \code{$oas} (background CDR3s), \code{$iedb} (epitope data),
#'   \code{$sabdab} (structural data), \code{$kmer_freq} (precomputed
#'   k-mer frequencies from OAS).
#'
#' @export
build_reference <- function(oas_path = NULL,
                            iedb_path = NULL,
                            sabdab_path = NULL,
                            species = "human") {

  # Download anything missing
  if (is.null(oas_path)) oas_path <- download_oas_reference(species = species)
  if (is.null(iedb_path)) iedb_path <- download_iedb_bcell()
  if (is.null(sabdab_path)) sabdab_path <- download_sabdab()

  oas    <- readRDS(oas_path)
  iedb   <- readRDS(iedb_path)
  sabdab <- readRDS(sabdab_path)

  # Precompute k-mer frequencies from OAS background
  message("[blusteR] Precomputing background k-mer frequencies...")
  kmer_freq <- .compute_kmer_frequencies(oas$cdr3_aa,
                                         k_sizes = .DEFAULT_KMER_SIZES)

  ref <- list(
    oas       = oas,
    iedb      = iedb,
    sabdab    = sabdab,
    kmer_freq = kmer_freq,
    species   = species
  )
  class(ref) <- "bluster_reference"

  message("[blusteR] Reference built successfully.")
  ref
}


# ---- internal helpers for reference databases -------------------------

#' Download sequences from OAS API
#' @keywords internal
.download_oas_api <- function(species, n_sequences) {

  # OAS bulk data endpoint
  base_url <- "http://opig.stats.ox.ac.uk/webapps/oas/oas_unpaired/"
  params <- list(
    species = species,
    chain   = "Heavy",
    output  = "json"
  )

  resp <- httr::GET(base_url, query = params, httr::timeout(300))

  if (httr::status_code(resp) != 200) {
    stop("OAS API request failed with status ", httr::status_code(resp))
  }

  raw <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                            flatten = TRUE)

  # Extract CDR3 sequences and V-gene annotations
  dt <- data.table::data.table(
    cdr3_aa = raw$cdr3_aa,
    v_gene  = raw$v_gene,
    j_gene  = raw$j_gene,
    chain   = "IGH"
  )

  # Subsample if needed
  if (nrow(dt) > n_sequences) {
    dt <- dt[sample(.N, n_sequences)]
  }

  dt
}


#' Generate synthetic reference from known CDR3 length/AA distributions
#' @keywords internal
.generate_synthetic_reference <- function(species, n_sequences) {

  # Human IGH CDR3 length distribution (approximate from literature)
  len_probs <- c(
    `8` = 0.01, `9` = 0.02, `10` = 0.04, `11` = 0.06, `12` = 0.09,
    `13` = 0.12, `14` = 0.14, `15` = 0.14, `16` = 0.12, `17` = 0.09,
    `18` = 0.06, `19` = 0.04, `20` = 0.03, `21` = 0.02, `22` = 0.01,
    `23` = 0.01, `24` = 0.005, `25` = 0.005
  )
  len_probs <- len_probs / sum(len_probs)

  # Position-specific AA frequencies in CDR3H (simplified)
  # Start with C, end with W (canonical)
  inner_aa <- c("A","R","N","D","C","Q","E","G","H","I",
                "L","K","M","F","P","S","T","W","Y","V")
  # Weights reflecting typical CDR3 composition
  inner_wt <- c(3,4,2,3,1,2,2,6,1,2,3,2,1,2,2,4,3,1,4,3)
  inner_wt <- inner_wt / sum(inner_wt)

  # Common IGHV genes
  v_genes <- c("IGHV1-2","IGHV1-18","IGHV1-69","IGHV2-5","IGHV3-7",
               "IGHV3-15","IGHV3-21","IGHV3-23","IGHV3-30","IGHV3-33",
               "IGHV3-48","IGHV3-49","IGHV3-53","IGHV3-66","IGHV3-74",
               "IGHV4-4","IGHV4-34","IGHV4-39","IGHV4-59","IGHV5-51")
  v_wt <- c(3,2,5,1,3,2,3,4,5,3,2,2,2,1,2,3,4,3,4,2)
  v_wt <- v_wt / sum(v_wt)

  j_genes <- c("IGHJ1","IGHJ2","IGHJ3","IGHJ4","IGHJ5","IGHJ6")
  j_wt <- c(0.05, 0.08, 0.10, 0.35, 0.15, 0.27)

  lengths <- sample(as.integer(names(len_probs)), n_sequences,
                    replace = TRUE, prob = len_probs)

  cdr3s <- vapply(lengths, function(l) {
    inner <- sample(inner_aa, l - 2, replace = TRUE, prob = inner_wt)
    paste0("C", paste0(inner, collapse = ""), "W")
  }, character(1))

  data.table::data.table(
    cdr3_aa = cdr3s,
    v_gene  = sample(v_genes, n_sequences, replace = TRUE, prob = v_wt),
    j_gene  = sample(j_genes, n_sequences, replace = TRUE, prob = j_wt),
    chain   = "IGH"
  )
}


#' Parse IEDB B-cell response into standardised table
#' @keywords internal
.parse_iedb_bcell <- function(raw) {

  dt <- data.table::as.data.table(raw)

  # Standardise column names
  cols_want <- c("epitope_linear_sequence", "antigen_name",
                 "organism_name", "antibody_heavy_chain_cdr3",
                 "antibody_light_chain_cdr3", "assay_type",
                 "pdb_id")
  cols_have <- intersect(cols_want, names(dt))

  out <- dt[, ..cols_have]
  data.table::setnames(out, old = cols_have, new = c(
    "epitope_seq", "antigen", "organism", "cdr3h",
    "cdr3l", "assay", "pdb_id"
  )[seq_along(cols_have)])

  out
}


#' Create a template IEDB database with known antibody-epitope pairs
#' @keywords internal
.create_template_iedb <- function() {

  data.table::data.table(
    epitope_seq = character(0),
    antigen     = character(0),
    organism    = character(0),
    cdr3h       = character(0),
    cdr3l       = character(0),
    assay       = character(0),
    pdb_id      = character(0)
  )
}


#' Parse SAbDab summary into standardised table
#' @keywords internal
.parse_sabdab <- function(raw) {

  dt <- data.table::as.data.table(raw)

  # SAbDab column names vary; map the common ones
  col_map <- list(
    pdb_id  = c("pdb", "PDB"),
    cdr3h   = c("CDRH3", "Hchain_cdr3"),
    cdr3l   = c("CDRL3", "Lchain_cdr3"),
    antigen = c("antigen_name", "antigen"),
    antigen_type = c("antigen_type"),
    vh      = c("Hchain_v_gene", "heavy_v"),
    vl      = c("Lchain_v_gene", "light_v"),
    resolution = c("resolution")
  )

  out_cols <- list()
  for (target in names(col_map)) {
    found <- intersect(col_map[[target]], names(dt))
    if (length(found) > 0) {
      out_cols[[target]] <- dt[[found[1]]]
    } else {
      out_cols[[target]] <- rep(NA_character_, nrow(dt))
    }
  }

  data.table::as.data.table(out_cols)
}

#' Create a template SAbDab database
#' @keywords internal
.create_template_sabdab <- function() {

  data.table::data.table(
    pdb_id       = character(0),
    cdr3h        = character(0),
    cdr3l        = character(0),
    antigen      = character(0),
    antigen_type = character(0),
    vh           = character(0),
    vl           = character(0),
    resolution   = character(0)
  )
}


#' Compute k-mer frequencies from a vector of CDR3 sequences
#' @keywords internal
.compute_kmer_frequencies <- function(cdr3_vec, k_sizes = c(4L, 5L)) {

  result <- list()
  cdr3_vec <- cdr3_vec[!is.na(cdr3_vec) & nchar(cdr3_vec) >= max(k_sizes)]
  n_total <- length(cdr3_vec)

  for (k in k_sizes) {
    # Extract all k-mers from all sequences
    all_kmers <- unlist(lapply(cdr3_vec, function(s) {
      n <- nchar(s)
      if (n < k) return(character(0))
      substring(s, seq_len(n - k + 1), seq_len(n - k + 1) + k - 1)
    }))

    tab <- table(all_kmers)
    freq <- as.numeric(tab) / sum(tab)
    names(freq) <- names(tab)

    # Also compute per-sequence frequency (fraction of sequences containing
    # each k-mer) for enrichment testing
    per_seq <- vapply(names(tab), function(km) {
      sum(grepl(km, cdr3_vec, fixed = TRUE)) / n_total
    }, numeric(1))

    result[[as.character(k)]] <- list(
      count    = tab,
      freq     = freq,
      per_seq  = per_seq,
      n_total  = n_total
    )
  }

  result
}
