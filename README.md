# blusteR: B-Cell Receptor Clustering by Shared Antigen Specificity

**Clustering BCR sequences by shared antigen-binding specificity — a
GLIPH2-inspired approach adapted for antibody biology.**

## Motivation

GLIPH2 (Huang et al., *Nature Biotechnology* 2020) revolutionised
TCR analysis by grouping receptors that recognise the same peptide-MHC
into "specificity groups."  No equivalent tool existed for B-cell
receptors, despite the growing volume of single-cell BCR-seq data
from 10X Genomics and other platforms.

blusteR fills this gap.  It discovers convergent CDR3 motifs and
near-identical CDR3 sequences across patients, then clusters them into
candidate specificity groups — sets of B cells that likely bind the
same antigenic epitope.

## Key adaptations from TCR → BCR

| Feature | GLIPH2 (TCR) | blusteR (BCR) |
|---------|-------------|-------------|
| CDR3 region | CDR3β (short, ~13 aa) | CDR3H (longer, ~15 aa) + optional CDR3L |
| K-mer size | 2–3 aa | 4–5 aa (captures longer motifs) |
| Distance metric | Hamming | BLOSUM62-weighted (SHM-aware) |
| Somatic hypermutation | N/A | Degenerate motifs + conservative-substitution tolerance |
| Chain pairing | CDR3α optional | Full heavy + light chain paired analysis |
| Reference background | Naive TCR repertoire | Observed Antibody Space (OAS) |
| MHC restriction | HLA enrichment test | N/A (antibodies bind antigen directly) |
| Clonal expansion | Not modelled | Clone collapsing to avoid inflation |
| Antigen databases | VDJdb (TCR epitopes) | IEDB B-cell epitopes + SAbDab structures |

## Reference databases

blusteR uses three public databases:

1. **Observed Antibody Space (OAS)** — Background CDR3 distribution
   for statistical enrichment testing.  OAS is the largest collection
   of antibody sequences from NGS studies (>2 billion sequences).
   *Source: opig.stats.ox.ac.uk/webapps/oas/*

2. **IEDB (Immune Epitope Database)** — Curated B-cell epitopes with
   antibody sequence annotations where available.  Used to annotate
   clusters with known antigen specificities.
   *Source: iedb.org*

3. **SAbDab (Structural Antibody Database)** — All antibody structures
   from the PDB with CDR3 and antigen annotations.  Provides
   structurally validated CDR3–antigen pairs for annotation.
   *Source: opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab*

## Installation

```r
# Install dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")

install.packages(c("data.table", "stringdist", "igraph",
                    "ggplot2", "jsonlite", "httr"))

# Install blusteR from local source
install.packages("/path/to/blusteR", repos = NULL, type = "source")

# Or with devtools from GitHub (once published)
# devtools::install_github("username/blusteR")
```

## Quick start

### From 10X Genomics VDJ output

```r
library(blusteR)

# Load BCR data from Cell Ranger VDJ
bcr <- load_10x_bcr("/path/to/cellranger_vdj/outs/")

# Run the full pipeline
result <- bluster(bcr)

# Explore results
summarize_clusters(result)
plot_cluster_network(result)
plot_cluster_summary(result)

# Export
export_clusters(result, prefix = "my_study", dir = "results/")
```

### From scRepertoire

```r
library(blusteR)
library(scRepertoire)

# Your standard scRepertoire workflow
contig_list <- list(read.csv("sample1_contigs.csv"),
                    read.csv("sample2_contigs.csv"))
combined <- combineBCR(contig_list, samples = c("Patient1", "Patient2"))

# Convert to blusteR format
bcr <- load_screpertoire(combined)

# Run with paired heavy+light chain analysis
result <- bluster(bcr, chain = "paired")
```

### From AIRR-format data

```r
bcr <- load_airr("rearrangements.tsv", sample_id = "my_sample")
result <- bluster(bcr)
```

### Multi-sample analysis

```r
# Load multiple donors/timepoints
bcr1 <- load_10x_bcr("donor1/", sample_id = "D1_pre")
bcr2 <- load_10x_bcr("donor2/", sample_id = "D2_pre")
bcr3 <- load_10x_bcr("donor1_post/", sample_id = "D1_post")

bcr_all <- rbind(bcr1, bcr2, bcr3)

# Cross-sample clusters provide strongest evidence for shared specificity
result <- bluster(bcr_all, chain = "heavy", annotate = TRUE)
```

## Algorithm overview

```
 ┌──────────────────────────────────────────────────────────────┐
 │                      INPUT                                   │
 │  10X VDJ / scRepertoire / AIRR  →  standardised data.table  │
 └───────────────────────┬──────────────────────────────────────┘
                         │
          ┌──────────────┴──────────────┐
          ▼                             ▼
 ┌─────────────────┐          ┌──────────────────┐
 │  LOCAL MOTIFS    │          │ GLOBAL SIMILARITY │
 │  (k=4,5 aa)     │          │  (BLOSUM62-       │
 │  Enrichment vs  │          │   weighted dist)  │
 │  OAS background │          │                   │
 └────────┬────────┘          └────────┬──────────┘
          │                            │
          └──────────┬─────────────────┘
                     ▼
          ┌─────────────────────┐
          │  NETWORK CLUSTERING │
          │  Connected          │
          │  components or      │
          │  Louvain detection  │
          └──────────┬──────────┘
                     ▼
          ┌─────────────────────┐
          │  CLUSTER SCORING    │
          │  • V-gene enrichment│
          │  • Clonal diversity │
          │  • Sample diversity │
          │  • CDR3 length SD   │
          │  • Composite score  │
          └──────────┬──────────┘
                     ▼
          ┌─────────────────────┐
          │ EPITOPE ANNOTATION  │
          │  IEDB + SAbDab      │
          │  CDR3 matching      │
          └──────────┬──────────┘
                     ▼
          ┌─────────────────────┐
          │       OUTPUT        │
          │  Clusters, motifs,  │
          │  network, plots     │
          └─────────────────────┘
```

## Cluster scoring

Each cluster receives a composite **blusteR score** (0–1) computed from:

| Component (weight) | Rationale |
|---------------------|-----------|
| Cluster size (20%) | Larger convergent groups are more statistically robust |
| Clonal diversity (25%) | High diversity = convergent selection, not clonal expansion |
| Sample diversity (25%) | Cross-patient convergence is strong evidence of shared specificity |
| V-gene enrichment (15%) | Shared V-gene usage supports structural similarity |
| CDR3 length homogeneity (15%) | Conserved length suggests conserved binding geometry |

## Output

`bluster()` returns a `bluster_result` object containing:

- `$clusters` — data.table of cluster-level summaries and scores
- `$membership` — data.table mapping each cell_id to its cluster
- `$edges` — combined edge table (local + global)
- `$graph` — igraph network object
- `$motifs` — enriched CDR3 motifs
- `$annotations` — IEDB/SAbDab epitope matches (if `annotate = TRUE`)

## Citation

If you use blusteR, please cite the tools it builds on:

- Olsen et al. "Observed Antibody Space: A diverse database..."
  *Protein Science* (2022)
- Vita et al. "The Immune Epitope Database (IEDB)..."
  *Nucleic Acids Research* (2019)
- Dunbar et al. "SAbDab: the structural antibody database."
  *Nucleic Acids Research* (2014)
- Huang et al. "Analyzing the Mycobacterium tuberculosis immune
  response by T-cell receptor clustering with GLIPH2..."
  *Nature Biotechnology* (2020)

## License

MIT
