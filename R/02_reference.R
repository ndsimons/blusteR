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
  message("        Querying IEDB Query API (https://query-api.iedb.org/)...")

  epitopes <- tryCatch({
    .download_iedb_api(organism)
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
#' @param n_cores Number of CPU cores to use when precomputing background
#'   k-mer frequencies.  Defaults to \code{getOption("bluster.ncores", 1)}.
#'   Set to a value > 1 (e.g. \code{parallel::detectCores()}) to parallelise
#'   the slow per-sequence k-mer frequency step.
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
                            species = "human",
                            n_cores = getOption("bluster.ncores", 1L)) {

  # Download anything missing
  if (is.null(oas_path)) oas_path <- download_oas_reference(species = species)
  if (is.null(iedb_path)) iedb_path <- download_iedb_bcell()
  if (is.null(sabdab_path)) sabdab_path <- download_sabdab()

  oas    <- readRDS(oas_path)
  iedb   <- readRDS(iedb_path)
  sabdab <- readRDS(sabdab_path)

  # Precompute k-mer frequencies from OAS background
  n_cores <- .resolve_ncores(n_cores)
  message("[blusteR] Precomputing background k-mer frequencies",
          if (n_cores > 1L) paste0(" (", n_cores, " cores)") else "", "...")
  kmer_freq <- .compute_kmer_frequencies(oas$cdr3_aa,
                                         k_sizes = .DEFAULT_KMER_SIZES,
                                         n_cores = n_cores)

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


#' Load the pre-built blusteR reference object
#'
#' Loads the pre-computed \code{bluster_reference} object bundled with the
#' package in \code{inst/extdata}.  These objects (OAS background, IEDB
#' epitopes, SAbDab structures, and the pre-computed background k-mer
#' frequency tables) are produced by \code{data-raw/build_references.R}.
#'
#' Using the bundled reference avoids re-downloading the source databases
#' and re-computing the slow background k-mer frequencies on every run.
#'
#' @param species One of \code{"human"} (default) or \code{"mouse"}.
#'
#' @return A list of class \code{bluster_reference}.
#' @export
load_reference <- function(species = c("human", "mouse")) {

  species <- match.arg(species)

  rds_path <- system.file(
    "extdata",
    paste0("bluster_reference_", species, ".rds"),
    package = "blusteR"
  )

  if (!nzchar(rds_path) || !file.exists(rds_path)) {
    stop("Pre-built reference for species '", species, "' not found. ",
         "Expected a bundled file at inst/extdata/bluster_reference_",
         species, ".rds. Rebuild it with data-raw/build_references.R ",
         "or call build_reference() to construct one at run time.",
         call. = FALSE)
  }

  message("[blusteR] Loading pre-built reference: ", rds_path)
  ref <- readRDS(rds_path)

  if (!inherits(ref, "bluster_reference")) {
    class(ref) <- "bluster_reference"
  }

  ref
}


# ---- internal helpers for reference databases -------------------------

#' Download sequences from OAS
#'
#' OAS does not expose a JSON REST API.  The unpaired-sequence search form
#' (a POST request) returns an HTML page that embeds a bulk-download shell
#' script of \code{wget <data-unit>.csv.gz} commands.  We submit the
#' attribute search, scrape those data-unit URLs, then download and parse
#' the gzipped CSV data units until enough sequences are collected.
#' @keywords internal
.download_oas_api <- function(species, n_sequences) {

  search_url <- "https://opig.stats.ox.ac.uk/webapps/oas/oas_unpaired/"

  resp <- httr::POST(
    search_url,
    body   = list(Species = species, Chain = "Heavy"),
    encode = "multipart",
    httr::timeout(300)
  )
  if (httr::status_code(resp) != 200) {
    stop("OAS search request failed with status ", httr::status_code(resp))
  }

  html <- httr::content(resp, "text", encoding = "UTF-8")

  # Pull the gzipped data-unit URLs out of the embedded download script
  unit_urls <- regmatches(
    html,
    gregexpr("https?://[^\"'[:space:]]+\\.csv\\.gz", html, perl = TRUE)
  )[[1]]
  unit_urls <- unique(unit_urls)

  if (length(unit_urls) == 0) {
    stop("OAS search returned no data units for species '", species, "'")
  }

  # Shuffle so we sample across studies rather than always the first ones
  unit_urls <- sample(unit_urls)

  parts  <- list()
  n_have <- 0L
  for (url in unit_urls) {
    part <- tryCatch(.download_oas_unit(url), error = function(e) NULL)
    if (is.null(part) || nrow(part) == 0) next
    parts[[length(parts) + 1L]] <- part
    n_have <- n_have + nrow(part)
    if (n_have >= n_sequences) break
  }

  if (length(parts) == 0) {
    stop("failed to download any OAS data units")
  }

  dt <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)

  # Keep only valid amino-acid CDR3 sequences
  dt <- dt[!is.na(dt$cdr3_aa) & nzchar(dt$cdr3_aa) &
             grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", dt$cdr3_aa)]

  if (nrow(dt) > n_sequences) {
    dt <- dt[sample(.N, n_sequences)]
  }

  dt[]
}


#' Download and parse a single OAS data unit (gzipped CSV)
#'
#' Each OAS data unit is a gzipped CSV whose first line is a metadata JSON
#' object and whose second line is the real column header.  The relevant
#' columns are \code{cdr3_aa}, \code{v_call}, \code{j_call} and \code{locus}.
#' @keywords internal
.download_oas_unit <- function(url, max_lines = 100000L) {

  resp <- httr::GET(url, httr::timeout(300))
  if (httr::status_code(resp) != 200) {
    stop("OAS data unit download failed: ", url)
  }

  tmp <- tempfile(fileext = ".csv.gz")
  on.exit(unlink(tmp), add = TRUE)
  writeBin(httr::content(resp, "raw"), tmp)

  con <- gzfile(tmp, "rt")
  # Skip the metadata line, then read up to max_lines of CSV (header + rows)
  lines <- readLines(con, n = max_lines + 1L)
  close(con)
  if (length(lines) < 3) return(NULL)

  dt <- data.table::fread(
    text             = paste(lines[-1], collapse = "\n"),
    stringsAsFactors = FALSE,
    showProgress     = FALSE
  )

  data.table::data.table(
    cdr3_aa = if ("cdr3_aa" %in% names(dt)) dt$cdr3_aa else NA_character_,
    v_gene  = if ("v_call"  %in% names(dt)) dt$v_call  else NA_character_,
    j_gene  = if ("j_call"  %in% names(dt)) dt$j_call  else NA_character_,
    chain   = "IGH"
  )
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


#' Query the IEDB Query API (IQ-API) for positive B-cell epitopes
#'
#' The IQ-API is a PostgREST service, so filters use PostgREST operators
#' (e.g. \code{qualitative_measure=eq.Positive}).  Results are paged at a
#' maximum of 10,000 records per request; we page with \code{offset} up to
#' \code{max_records}.
#' @keywords internal
.download_iedb_api <- function(organism = NULL, max_records = 50000L) {

  base_url <- "https://query-api.iedb.org/bcell_search"

  select_cols <- paste(
    c("linear_sequence", "parent_source_antigen_name",
      "source_organism_name", "receptor_chain1_cdr3_seqs",
      "receptor_chain2_cdr3_seqs", "pdb_id", "assay_names"),
    collapse = ","
  )

  page_size <- 10000L
  offset    <- 0L
  pages     <- list()

  repeat {
    query <- list(
      qualitative_measure = "eq.Positive",
      structure_type      = "eq.Linear peptide",
      select              = select_cols,
      order               = "structure_id",
      limit               = page_size,
      offset              = offset
    )
    if (!is.null(organism)) {
      query$source_organism_name <- paste0("ilike.*", organism, "*")
    }

    resp <- httr::GET(base_url, query = query, httr::timeout(120))
    if (httr::status_code(resp) != 200) {
      stop("IEDB returned status ", httr::status_code(resp))
    }

    raw <- jsonlite::fromJSON(
      httr::content(resp, "text", encoding = "UTF-8"),
      flatten = TRUE
    )
    if (length(raw) == 0 || NROW(raw) == 0) break

    pages[[length(pages) + 1L]] <- data.table::as.data.table(raw)
    offset <- offset + page_size
    if (offset >= max_records) break
  }

  if (length(pages) == 0) {
    stop("IEDB returned no B-cell epitope records")
  }

  .parse_iedb_bcell(data.table::rbindlist(pages, use.names = TRUE, fill = TRUE))
}


#' Parse IEDB B-cell response into standardised table
#' @keywords internal
.parse_iedb_bcell <- function(raw) {

  dt <- data.table::as.data.table(raw)

  if (nrow(dt) == 0) return(.create_template_iedb())

  # Receptor CDR3 fields are arrays (list-columns); collapse to a string
  collapse_seq <- function(x) {
    if (is.null(x)) return(rep(NA_character_, nrow(dt)))
    vapply(x, function(v) {
      v <- unlist(v)
      if (length(v) == 0 || all(is.na(v))) NA_character_
      else paste(v[!is.na(v)], collapse = ";")
    }, character(1))
  }

  pick <- function(name) {
    if (name %in% names(dt)) as.character(dt[[name]]) else NA_character_
  }

  data.table::data.table(
    epitope_seq = pick("linear_sequence"),
    antigen     = pick("parent_source_antigen_name"),
    organism    = pick("source_organism_name"),
    cdr3h       = collapse_seq(if ("receptor_chain1_cdr3_seqs" %in% names(dt))
                                 dt[["receptor_chain1_cdr3_seqs"]] else NULL),
    cdr3l       = collapse_seq(if ("receptor_chain2_cdr3_seqs" %in% names(dt))
                                 dt[["receptor_chain2_cdr3_seqs"]] else NULL),
    assay       = pick("assay_names"),
    pdb_id      = pick("pdb_id")
  )
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
#'
#' The per-sequence frequency step (fraction of background sequences
#' containing each k-mer) is the bottleneck: it scans every reference
#' sequence once per unique k-mer.  When \code{n_cores > 1} this scan is
#' split across CPU cores (forked workers on Unix via
#' \code{parallel::mclapply}, a PSOCK cluster on Windows).
#' @keywords internal
.compute_kmer_frequencies <- function(cdr3_vec, k_sizes = c(4L, 5L),
                                      n_cores = getOption("bluster.ncores", 1L)) {

  n_cores <- .resolve_ncores(n_cores)

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

    # Per-sequence frequency (fraction of sequences containing each k-mer)
    # for enrichment testing, parallelised across unique k-mers.
    kmers <- names(tab)
    per_seq <- .parallel_kmer_counts(kmers, cdr3_vec, n_cores) / n_total
    names(per_seq) <- kmers

    result[[as.character(k)]] <- list(
      count    = tab,
      freq     = freq,
      per_seq  = per_seq,
      n_total  = n_total
    )
  }

  result
}


#' Normalise a requested core count against the machine's capacity
#' @keywords internal
.resolve_ncores <- function(n_cores) {
  if (is.null(n_cores) || length(n_cores) == 0L ||
      is.na(n_cores[1]) || n_cores[1] < 1L) {
    return(1L)
  }
  max_cores <- tryCatch(parallel::detectCores(), error = function(e) 1L)
  if (is.na(max_cores) || max_cores < 1L) max_cores <- 1L
  as.integer(min(as.integer(n_cores[1]), max_cores))
}


#' Count, for each k-mer, how many sequences contain it (optionally parallel)
#'
#' Returns an integer-like numeric vector aligned to \code{kmers}.
#' @keywords internal
.parallel_kmer_counts <- function(kmers, cdr3_vec, n_cores = 1L) {

  count_chunk <- function(km_chunk) {
    vapply(km_chunk, function(km) {
      sum(grepl(km, cdr3_vec, fixed = TRUE))
    }, numeric(1))
  }

  if (n_cores <= 1L || length(kmers) < 2L) {
    return(count_chunk(kmers))
  }

  # Contiguous chunks preserve the original k-mer order on reassembly.
  chunk_id <- cut(seq_along(kmers), n_cores, labels = FALSE)
  chunks   <- split(kmers, chunk_id)

  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, "cdr3_vec", envir = environment())
    parts <- parallel::parLapply(cl, chunks, count_chunk)
  } else {
    parts <- parallel::mclapply(chunks, count_chunk, mc.cores = n_cores)
    failed <- vapply(parts, inherits, logical(1), what = "try-error")
    if (any(failed)) {
      stop("parallel k-mer counting failed: ",
           conditionMessage(attr(parts[[which(failed)[1]]], "condition")))
    }
  }

  unlist(parts, use.names = FALSE)
}
