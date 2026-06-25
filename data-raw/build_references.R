#!/usr/bin/env Rscript

# =====================================================================
# build_references.R  --  STAND-ALONE
# ---------------------------------------------------------------------
# Self-contained script to download the OAS / IEDB / SAbDab source data
# and build the pre-computed blusteR reference objects (including the
# slow background k-mer frequency tables).
#
# This script has NO dependency on the blusteR package being installed
# or present. Copy it anywhere (e.g. onto a VM) and run it. It only
# needs these CRAN packages: data.table, httr, jsonlite (parallel and
# stats ship with R).
#
# It writes, per species:
#       <out>/bluster_reference_<species>.rds
# where <out> defaults to the current directory. Copy the resulting
# .rds files into the package's inst/extdata/ to have them bundled and
# loaded automatically at run time.
#
# Usage:
#   Rscript build_references.R
#   Rscript build_references.R --species=human
#   Rscript build_references.R --species=human,mouse --n=500000
#   Rscript build_references.R --cores=8 --out=./refs --force
#
# Options:
#   --species=...  Comma-separated species to build (default: human,mouse)
#   --n=...        Max OAS background sequences per species (default: 500000)
#   --cores=...    Cores for k-mer precomputation (default: detected - 1)
#   --out=...      Output directory for the .rds files (default: .)
#   --cache=...    Directory for intermediate downloads (default: <out>/cache)
#   --force        Re-download source data and rebuild even if cached
# =====================================================================

# ---- Dependencies ----------------------------------------------------

for (dep in c("data.table", "httr", "jsonlite")) {
  if (!requireNamespace(dep, quietly = TRUE)) {
    stop("Required package '", dep, "' is not installed. Install with:\n",
         "  install.packages(c(\"data.table\", \"httr\", \"jsonlite\"))",
         call. = FALSE)
  }
}

DEFAULT_KMER_SIZES <- c(4L, 5L)

# =====================================================================
# Reference-building functions (embedded; no package required)
# =====================================================================

# ---- OAS background --------------------------------------------------

download_oas_reference <- function(dest_dir,
                                   species     = "human",
                                   n_sequences = 500000L,
                                   force       = FALSE) {

  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  rds_path <- file.path(dest_dir, paste0("oas_ref_", species, ".rds"))

  if (file.exists(rds_path) && !force) {
    message("[blusteR] OAS reference already exists: ", rds_path)
    return(invisible(rds_path))
  }

  message("[blusteR] Building OAS reference for ", species, "...")
  message("        Querying OAS (https://opig.stats.ox.ac.uk/webapps/oas/)...")

  ref <- tryCatch({
    .download_oas_api(species, n_sequences)
  }, error = function(e) {
    message("[blusteR] OAS API unavailable: ", conditionMessage(e))
    message("        Generating synthetic reference from known distributions...")
    .generate_synthetic_reference(species, n_sequences)
  })

  saveRDS(ref, rds_path)
  message("[blusteR] OAS reference saved: ", rds_path,
          " (", nrow(ref), " sequences)")
  invisible(rds_path)
}

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

  unit_urls <- regmatches(
    html,
    gregexpr("https?://[^\"'[:space:]]+\\.csv\\.gz", html, perl = TRUE)
  )[[1]]
  unit_urls <- unique(unit_urls)

  if (length(unit_urls) == 0) {
    stop("OAS search returned no data units for species '", species, "'")
  }

  unit_urls <- sample(unit_urls)

  parts  <- list()
  n_have <- 0L
  for (url in unit_urls) {
    part <- tryCatch(.download_oas_unit(url), error = function(e) NULL)
    if (is.null(part) || nrow(part) == 0) next
    parts[[length(parts) + 1L]] <- part
    n_have <- n_have + nrow(part)
    message("        ... ", n_have, " sequences collected")
    if (n_have >= n_sequences) break
  }

  if (length(parts) == 0) stop("failed to download any OAS data units")

  dt <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  dt <- dt[!is.na(dt$cdr3_aa) & nzchar(dt$cdr3_aa) &
             grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", dt$cdr3_aa)]

  if (nrow(dt) > n_sequences) dt <- dt[sample(.N, n_sequences)]
  dt[]
}

.download_oas_unit <- function(url, max_lines = 100000L) {

  resp <- httr::GET(url, httr::timeout(300))
  if (httr::status_code(resp) != 200) {
    stop("OAS data unit download failed: ", url)
  }

  tmp <- tempfile(fileext = ".csv.gz")
  on.exit(unlink(tmp), add = TRUE)
  writeBin(httr::content(resp, "raw"), tmp)

  con <- gzfile(tmp, "rt")
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

.generate_synthetic_reference <- function(species, n_sequences) {

  len_probs <- c(
    `8` = 0.01, `9` = 0.02, `10` = 0.04, `11` = 0.06, `12` = 0.09,
    `13` = 0.12, `14` = 0.14, `15` = 0.14, `16` = 0.12, `17` = 0.09,
    `18` = 0.06, `19` = 0.04, `20` = 0.03, `21` = 0.02, `22` = 0.01,
    `23` = 0.01, `24` = 0.005, `25` = 0.005
  )
  len_probs <- len_probs / sum(len_probs)

  inner_aa <- c("A","R","N","D","C","Q","E","G","H","I",
                "L","K","M","F","P","S","T","W","Y","V")
  inner_wt <- c(3,4,2,3,1,2,2,6,1,2,3,2,1,2,2,4,3,1,4,3)
  inner_wt <- inner_wt / sum(inner_wt)

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

# ---- IEDB B-cell epitopes -------------------------------------------

download_iedb_bcell <- function(dest_dir, organism = NULL, force = FALSE) {

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

  if (length(pages) == 0) stop("IEDB returned no B-cell epitope records")

  .parse_iedb_bcell(data.table::rbindlist(pages, use.names = TRUE, fill = TRUE))
}

.parse_iedb_bcell <- function(raw) {

  dt <- data.table::as.data.table(raw)
  if (nrow(dt) == 0) return(.create_template_iedb())

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

.create_template_iedb <- function() {
  data.table::data.table(
    epitope_seq = character(0), antigen = character(0),
    organism    = character(0), cdr3h   = character(0),
    cdr3l       = character(0), assay   = character(0),
    pdb_id      = character(0)
  )
}

# ---- SAbDab structures ----------------------------------------------

download_sabdab <- function(dest_dir, force = FALSE) {

  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  rds_path <- file.path(dest_dir, "sabdab.rds")

  if (file.exists(rds_path) && !force) {
    message("[blusteR] SAbDab data already exists: ", rds_path)
    return(invisible(rds_path))
  }

  message("[blusteR] Downloading SAbDab structural antibody data...")

  sabdab <- tryCatch({
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

.parse_sabdab <- function(raw) {

  dt <- data.table::as.data.table(raw)

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
    out_cols[[target]] <- if (length(found) > 0) dt[[found[1]]]
                          else rep(NA_character_, nrow(dt))
  }

  data.table::as.data.table(out_cols)
}

.create_template_sabdab <- function() {
  data.table::data.table(
    pdb_id       = character(0), cdr3h        = character(0),
    cdr3l        = character(0), antigen      = character(0),
    antigen_type = character(0), vh           = character(0),
    vl           = character(0), resolution   = character(0)
  )
}

# ---- k-mer frequency precomputation (parallel) ----------------------

.resolve_ncores <- function(n_cores) {
  if (is.null(n_cores) || length(n_cores) == 0L ||
      is.na(n_cores[1]) || n_cores[1] < 1L) {
    return(1L)
  }
  max_cores <- tryCatch(parallel::detectCores(), error = function(e) 1L)
  if (is.na(max_cores) || max_cores < 1L) max_cores <- 1L
  as.integer(min(as.integer(n_cores[1]), max_cores))
}

.parallel_kmer_counts <- function(kmers, cdr3_vec, n_cores = 1L) {

  count_chunk <- function(km_chunk) {
    vapply(km_chunk, function(km) {
      sum(grepl(km, cdr3_vec, fixed = TRUE))
    }, numeric(1))
  }

  if (n_cores <= 1L || length(kmers) < 2L) {
    return(count_chunk(kmers))
  }

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

.compute_kmer_frequencies <- function(cdr3_vec, k_sizes = c(4L, 5L),
                                      n_cores = 1L) {

  n_cores <- .resolve_ncores(n_cores)

  result <- list()
  cdr3_vec <- cdr3_vec[!is.na(cdr3_vec) & nchar(cdr3_vec) >= max(k_sizes)]
  n_total <- length(cdr3_vec)

  for (k in k_sizes) {
    all_kmers <- unlist(lapply(cdr3_vec, function(s) {
      n <- nchar(s)
      if (n < k) return(character(0))
      substring(s, seq_len(n - k + 1), seq_len(n - k + 1) + k - 1)
    }))

    tab <- table(all_kmers)
    freq <- as.numeric(tab) / sum(tab)
    names(freq) <- names(tab)

    kmers <- names(tab)
    per_seq <- .parallel_kmer_counts(kmers, cdr3_vec, n_cores) / n_total
    names(per_seq) <- kmers

    result[[as.character(k)]] <- list(
      count = tab, freq = freq, per_seq = per_seq, n_total = n_total
    )
  }

  result
}

# ---- assemble the reference object ----------------------------------

build_reference <- function(oas_path, iedb_path, sabdab_path,
                            species = "human", n_cores = 1L) {

  oas    <- readRDS(oas_path)
  iedb   <- readRDS(iedb_path)
  sabdab <- readRDS(sabdab_path)

  n_cores <- .resolve_ncores(n_cores)
  message("[blusteR] Precomputing background k-mer frequencies",
          if (n_cores > 1L) paste0(" (", n_cores, " cores)") else "", "...")
  kmer_freq <- .compute_kmer_frequencies(oas$cdr3_aa,
                                         k_sizes = DEFAULT_KMER_SIZES,
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

# =====================================================================
# Command-line driver
# =====================================================================

args <- commandArgs(trailingOnly = TRUE)

get_opt <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
has_flag <- function(name) any(args == paste0("--", name))

species_arg  <- get_opt("species", "human,mouse")
species_list <- trimws(strsplit(species_arg, ",")[[1]])
species_list <- species_list[nzchar(species_list)]

n_sequences <- as.integer(get_opt("n", "500000"))
force       <- has_flag("force")

out_dir   <- normalizePath(get_opt("out", "."), mustWork = FALSE)
cache_dir <- normalizePath(get_opt("cache", file.path(out_dir, "cache")),
                           mustWork = FALSE)

default_cores <- max(1L, parallel::detectCores() - 1L)
n_cores <- as.integer(get_opt("cores", as.character(default_cores)))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L

dir.create(out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

message("[build_references] Species      : ", paste(species_list, collapse = ", "))
message("[build_references] Max sequences: ", n_sequences)
message("[build_references] Cores        : ", n_cores)
message("[build_references] Force rebuild: ", force)
message("[build_references] Cache dir    : ", cache_dir)
message("[build_references] Output dir   : ", out_dir)

# IEDB and SAbDab are species-independent; download once and reuse.
message("\n[build_references] === Downloading IEDB B-cell epitopes ===")
iedb_path <- download_iedb_bcell(dest_dir = cache_dir, force = force)

message("\n[build_references] === Downloading SAbDab structures ===")
sabdab_path <- download_sabdab(dest_dir = cache_dir, force = force)

build_one <- function(species) {
  message("\n[build_references] ======================================")
  message("[build_references] Building reference for: ", species)
  message("[build_references] ======================================")

  message("\n[build_references] --- Downloading OAS background (", species, ") ---")
  oas_path <- download_oas_reference(
    dest_dir    = cache_dir,
    species     = species,
    n_sequences = n_sequences,
    force       = force
  )

  message("\n[build_references] --- Assembling reference + k-mer frequencies ---")
  ref <- build_reference(
    oas_path    = oas_path,
    iedb_path   = iedb_path,
    sabdab_path = sabdab_path,
    species     = species,
    n_cores     = n_cores
  )

  out_path <- file.path(out_dir, paste0("bluster_reference_", species, ".rds"))
  saveRDS(ref, out_path, compress = "xz")

  size_mb <- round(file.info(out_path)$size / 1024^2, 2)
  message("[build_references] Saved ", out_path, " (", size_mb, " MB)")
  out_path
}

saved <- vapply(species_list, build_one, character(1))

message("\n[build_references] === Done ===")
for (p in saved) message("  - ", p)
message("\nCopy these .rds files into the package's inst/extdata/ so they are")
message("bundled and loaded automatically by build_reference() at run time.")
