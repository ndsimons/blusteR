library(testthat)
library(data.table)

# ---- helper: create mock BCR data ------------------------------------

make_mock_bcr <- function(n_cells = 50, n_samples = 2) {

  aa <- c("A","R","N","D","C","Q","E","G","H","I",
          "L","K","M","F","P","S","T","W","Y","V")

  v_genes <- c("IGHV3-30", "IGHV3-23", "IGHV4-34", "IGHV1-69", "IGHV3-7")
  j_genes <- c("IGHJ4", "IGHJ6", "IGHJ5")
  vl_genes <- c("IGKV1-39", "IGKV3-20", "IGLV1-51")
  jl_genes <- c("IGKJ1", "IGKJ2", "IGLJ2")
  samples  <- paste0("S", seq_len(n_samples))

  # Create some convergent CDR3s to ensure clusters form
  shared_core <- c("CARDGYSSGWY", "CARGYSSGWY", "CARDGYSSGWF",  # cluster A
                    "CAREVVAAGTDY", "CAREVVAAGTEY")               # cluster B

  rows <- list()
  for (i in seq_len(n_cells)) {
    sid <- sample(samples, 1)
    cid <- paste0(sid, "_cell", i)

    # Heavy chain
    if (i <= length(shared_core)) {
      cdr3h <- paste0(shared_core[i], sample(aa, 3, replace = TRUE) |> paste0(collapse = ""))
    } else {
      cdr3h <- paste0("C", paste0(sample(aa, sample(10:18, 1), replace = TRUE), collapse = ""), "W")
    }

    rows[[length(rows) + 1]] <- data.table(
      cell_id    = cid,
      sample_id  = sid,
      chain      = "IGH",
      chain_type = "heavy",
      v_gene     = sample(v_genes, 1),
      d_gene     = NA_character_,
      j_gene     = sample(j_genes, 1),
      c_gene     = NA_character_,
      cdr3_aa    = cdr3h,
      cdr3_nt    = NA_character_,
      clone_id   = paste0("clone_", ceiling(i / 3)),
      is_cell    = TRUE,
      productive = TRUE
    )

    # Light chain
    cdr3l <- paste0("C", paste0(sample(aa, sample(8:12, 1), replace = TRUE), collapse = ""), "F")
    chain_name <- sample(c("IGK", "IGL"), 1)
    rows[[length(rows) + 1]] <- data.table(
      cell_id    = cid,
      sample_id  = sid,
      chain      = chain_name,
      chain_type = "light",
      v_gene     = sample(vl_genes, 1),
      d_gene     = NA_character_,
      j_gene     = sample(jl_genes, 1),
      c_gene     = NA_character_,
      cdr3_aa    = cdr3l,
      cdr3_nt    = NA_character_,
      clone_id   = paste0("clone_", ceiling(i / 3)),
      is_cell    = TRUE,
      productive = TRUE
    )
  }

  rbindlist(rows, fill = TRUE)
}

# ---- tests -----------------------------------------------------------

test_that("mock data has correct structure", {
  bcr <- make_mock_bcr()
  expect_true(is.data.table(bcr))
  expect_true(all(c("cell_id", "chain_type", "cdr3_aa") %in% names(bcr)))
  expect_gt(nrow(bcr), 0)
  expect_true("heavy" %in% bcr$chain_type)
  expect_true("light" %in% bcr$chain_type)
})

test_that("filter_bcr removes short and ambiguous CDR3s", {
  bcr <- make_mock_bcr()

  # Add a bad row
  bad_row <- bcr[1]
  bad_row$cdr3_aa <- "CX"
  bad_row$cell_id <- "bad_cell"
  bcr <- rbind(bcr, bad_row)

  filtered <- filter_bcr(bcr, min_cdr3_len = 5)
  expect_lt(nrow(filtered), nrow(bcr))
  expect_false(any(grepl("X", filtered$cdr3_aa)))
})

test_that("BLOSUM62 matrix loads correctly", {
  b <- .get_blosum62()
  expect_true(is.matrix(b))
  expect_equal(nrow(b), 20)
  expect_equal(ncol(b), 20)
  expect_equal(b["A", "A"], 4)
  expect_true(isSymmetric(b))
})

test_that(".blosum_distance returns 0 for identical sequences", {
  b <- .get_blosum62()
  s <- strsplit("CARDGY", "")[[1]]
  expect_equal(.blosum_distance(s, s, b), 0)
})

test_that(".blosum_distance penalises non-conservative substitutions", {
  b <- .get_blosum62()
  s1 <- strsplit("CARDGY", "")[[1]]  # original
  s2 <- strsplit("CARDSY", "")[[1]]  # G→S: BLOSUM62 = 0 → 0.5 penalty
  s3 <- strsplit("CARDWY", "")[[1]]  # G→W: BLOSUM62 = -2 → 1 penalty

  d_conservative <- .blosum_distance(s1, s2, b)
  d_radical      <- .blosum_distance(s1, s3, b)

  expect_lt(d_conservative, d_radical)
})

test_that(".extract_kmers_per_seq works", {
  seqs <- c("ABCDEF", "GHIJKL")
  kmers <- .extract_kmers_per_seq(seqs, 4)
  expect_length(kmers, 2)
  expect_true("ABCD" %in% kmers[[1]])
  expect_true("BCDE" %in% kmers[[1]])
})

test_that(".collapse_clones reduces cell count", {
  bcr <- make_mock_bcr(n_cells = 30)
  collapsed <- .collapse_clones(bcr)
  expect_lte(nrow(collapsed), nrow(bcr))
})

test_that(".standardise_chain maps correctly", {
  expect_equal(.standardise_chain(c("IGH", "Heavy", "IGK", "Lambda")),
               c("IGH", "IGH", "IGK", "IGL"))
})

test_that("synthetic reference is generated correctly", {
  ref <- .generate_synthetic_reference("human", 1000)
  expect_true(is.data.table(ref))
  expect_equal(nrow(ref), 1000)
  expect_true(all(grepl("^C.*W$", ref$cdr3_aa)))  # canonical CDR3H
  expect_true(all(grepl("^IGHV", ref$v_gene)))
})

test_that("kmer frequency computation works", {
  seqs <- c("CARDGYSSGWY", "CAREVVAAGTDY", "CARDGYSSGWF", "CXYZABCDEFGW")
  freqs <- .compute_kmer_frequencies(seqs, k_sizes = c(4L))

  expect_true("4" %in% names(freqs))
  expect_true("CARD" %in% names(freqs[["4"]]$count))
  expect_gt(freqs[["4"]]$n_total, 0)
})

test_that(".compute_bluster_score returns values in [0, 1]", {
  scores <- .compute_bluster_score(
    n_members        = c(3, 10, 50),
    clonal_diversity = c(1.0, 0.8, 0.5),
    sample_diversity = c(0.5, 1.0, 0.2),
    vh_pvalue        = c(0.01, 0.001, 0.5),
    len_sd           = c(0.5, 1.0, 3.0)
  )
  expect_true(all(scores >= 0 & scores <= 1))
})

test_that("empty edge table has correct structure", {
  e <- .empty_edges()
  expect_true(is.data.table(e))
  expect_equal(nrow(e), 0)
  expect_true(all(c("seq_id_1", "seq_id_2", "distance") %in% names(e)))
})
