#' blusteR: B-Cell Lymphocyte Interaction by Paratope Hotspots
#'
#' Clusters BCR sequences by shared antigen-binding specificity.
#' Adapts the GLIPH2 paradigm to antibody biology, accounting for
#' somatic hypermutation, heavy/light chain pairing, and longer CDR3
#' regions.
#'
#' @section Key differences from GLIPH2 (TCR):
#' \itemize{
#'   \item Uses longer k-mers (4-5aa) for motif discovery due to longer CDR3H
#'   \item SHM-aware distance scoring (BLOSUM62-weighted global similarity)
#'   \item Paired heavy + light chain analysis when data is available
#'   \item Reference: Observed Antibody Space (OAS) instead of TCR reference
#'   \item Antigen annotation via IEDB B-cell epitopes and SAbDab structures
#'   \item No HLA enrichment (antibodies bind antigen directly)
#'   \item Clonal lineage awareness to avoid inflating cluster membership
#' }
#'
#' @docType package
#' @name blusteR-package
"_PACKAGE"

# ---- internal constants -----------------------------------------------

.blusteR_ENV <- new.env(parent = emptyenv())

# BLOSUM62 substitution scores for SHM-aware comparison
.BLOSUM62 <- NULL  # lazily loaded

# Amino acid properties for CDR3 physico-chemical grouping
.AA_GROUPS <- list(
  hydrophobic = c("A", "V", "I", "L", "M", "F", "W", "P"),
  polar       = c("S", "T", "N", "Q", "Y", "C"),
  positive    = c("R", "H", "K"),
  negative    = c("D", "E"),
  special     = c("G")
)

# Default k-mer sizes for BCR motif discovery (longer than TCR k=3)
.DEFAULT_KMER_SIZES <- c(4L, 5L)

# Default distance thresholds
.DEFAULT_GLOBAL_DIST <- 2L
.DEFAULT_LOCAL_MIN_FREQ <- 3L
.DEFAULT_PVALUE <- 0.05

#' Lazily load BLOSUM62 matrix
#' @keywords internal
.get_blosum62 <- function() {

  if (is.null(.BLOSUM62)) {
    aa <- c("A","R","N","D","C","Q","E","G","H","I",
            "L","K","M","F","P","S","T","W","Y","V")
    # Standard BLOSUM62 upper triangle, symmetric
    raw <- c(
       4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,
      -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,
      -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,
      -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,
       0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,
      -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,
      -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,
       0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,
      -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,
      -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,
      -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,
      -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,
      -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,
      -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,
      -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,
       1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,
       0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,
      -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,
      -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,
       0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4
    )
    m <- matrix(raw, nrow = 20, ncol = 20, dimnames = list(aa, aa))
    assign(".BLOSUM62", m, envir = parent.env(environment()))
  }
  .BLOSUM62
}
