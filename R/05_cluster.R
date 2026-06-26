#' Form and score blusteR clusters
#'
#' Combines local motif edges and global similarity edges into a
#' network, extracts connected components (or uses Louvain community
#' detection for large networks), then scores each cluster by:
#'
#' \enumerate{
#'   \item \strong{V-gene enrichment}: Fisher's exact test for over-
#'     representation of specific V genes within the cluster.
#'   \item \strong{CDR3 length homogeneity}: Whether cluster members
#'     share similar CDR3 lengths (suggesting structural convergence).
#'   \item \strong{Clonal diversity}: Number of distinct clonotypes
#'     contributing to the cluster (penalises clusters driven by
#'     clonal expansion rather than convergent selection).
#'   \item \strong{Sample diversity}: Number of distinct biological
#'     samples contributing (cross-sample convergence is stronger
#'     evidence of shared specificity).
#'   \item \strong{Composite blusteR score}: Weighted combination.
#' }
#'
#' @param bcr_data blusteR standard format data.table.
#' @param motif_results Output of \code{bluster_motifs()}.
#' @param global_edges Output of \code{bluster_global()}.
#' @param min_cluster_size Minimum members for a cluster to be retained
#'   (default 3).
#' @param clustering_method \code{"components"} (default) or
#'   \code{"louvain"} for community detection.
#' @param collapse_clones Logical; collapse cells from the same clone
#'   before clustering to avoid inflation (default TRUE).
#' @param n_cores Number of CPU cores for expanding motif edges (default
#'   \code{getOption("bluster.ncores", 1)}).  Values > 1 use forked
#'   workers; degrades to serial on Windows.
#'
#' @return A list of class \code{bluster_result}:
#' \describe{
#'   \item{clusters}{data.table with cluster-level summaries}
#'   \item{membership}{data.table mapping cell_id → cluster_id}
#'   \item{edges}{combined edge table}
#'   \item{graph}{igraph object}
#'   \item{motifs}{input motif results}
#'   \item{params}{parameters used}
#' }
#'
#' @export
bluster_cluster <- function(bcr_data,
                          motif_results,
                          global_edges,
                          min_cluster_size = 3L,
                          clustering_method = c("components", "louvain"),
                          collapse_clones = TRUE,
                          n_cores = getOption("bluster.ncores", 1L)) {

  clustering_method <- match.arg(clustering_method)
  n_cores <- .resolve_ncores(n_cores)

  # --- 1. Collapse clonally related cells --------------------------------
  if (collapse_clones) {
    bcr_data <- .collapse_clones(bcr_data)
  }

  # --- 2. Build combined edge list from motifs and global similarity -----
  edges <- .build_combined_edges(motif_results, global_edges, bcr_data, n_cores)

  if (nrow(edges) == 0) {
    message("[blusteR] No edges found; cannot form clusters.")
    return(.empty_bluster_result(bcr_data, motif_results))
  }

  # --- 3. Build graph and detect clusters --------------------------------
  g <- igraph::graph_from_data_frame(
    edges[, .(seq_id_1, seq_id_2)],
    directed = FALSE
  )

  if (clustering_method == "louvain" && igraph::vcount(g) > 50) {
    communities <- igraph::cluster_louvain(g)
    cluster_ids <- igraph::membership(communities)
  } else {
    comp <- igraph::components(g)
    cluster_ids <- comp$membership
  }

  # Assign cluster IDs
  membership <- data.table::data.table(
    cell_id    = names(cluster_ids),
    cluster_id = paste0("blusteR_", cluster_ids)
  )

  # --- 4. Filter by minimum size -----------------------------------------
  cluster_sizes <- membership[, .N, by = cluster_id]
  keep <- cluster_sizes[N >= min_cluster_size, cluster_id]
  membership <- membership[cluster_id %in% keep]

  if (nrow(membership) == 0) {
    message("[blusteR] No clusters meet minimum size threshold.")
    return(.empty_bluster_result(bcr_data, motif_results))
  }

  # --- 5. Score each cluster ---------------------------------------------
  message(sprintf("[blusteR] Scoring %d clusters...", length(keep)))

  cluster_info <- merge(membership, bcr_data, by = "cell_id", all.x = TRUE)

  clusters <- cluster_info[, {

    heavy <- .SD[chain_type == "heavy"]
    light <- .SD[chain_type == "light"]

    # V-gene enrichment (heavy chain)
    vgene_enrich <- .test_vgene_enrichment(
      heavy$v_gene, bcr_data[chain_type == "heavy"]$v_gene
    )

    # CDR3 length homogeneity
    len_sd <- if (nrow(heavy) > 1) sd(nchar(heavy$cdr3_aa), na.rm = TRUE) else 0

    # Clonal diversity
    n_clones <- length(unique(na.omit(heavy$clone_id)))

    # Sample diversity
    n_samples <- length(unique(.SD$sample_id))

    list(
      n_members          = .N,
      n_heavy            = nrow(heavy),
      n_light            = nrow(light),
      n_clones           = n_clones,
      n_samples          = n_samples,
      cdr3h_consensus    = .consensus_cdr3(heavy$cdr3_aa),
      cdr3h_len_mean     = mean(nchar(heavy$cdr3_aa), na.rm = TRUE),
      cdr3h_len_sd       = len_sd,
      top_vh             = .top_gene(heavy$v_gene),
      top_jh             = .top_gene(heavy$j_gene),
      vh_enrichment_p    = vgene_enrich$pvalue,
      vh_enriched_gene   = vgene_enrich$gene,
      clonal_diversity   = n_clones / max(.N, 1),
      sample_diversity   = n_samples / length(unique(bcr_data$sample_id))
    )
  }, by = cluster_id]

  # Compute composite blusteR score
  clusters[, bluster_score := .compute_bluster_score(
    n_members, clonal_diversity, sample_diversity,
    vh_enrichment_p, cdr3h_len_sd
  )]

  data.table::setorder(clusters, -bluster_score)

  message(sprintf("[blusteR] Clustering complete: %d clusters, %d cells assigned.",
                  nrow(clusters), nrow(membership)))

  result <- list(
    clusters   = clusters,
    membership = membership,
    edges      = edges,
    graph      = g,
    motifs     = motif_results,
    params     = list(
      min_cluster_size  = min_cluster_size,
      clustering_method = clustering_method,
      collapse_clones   = collapse_clones
    )
  )
  class(result) <- "bluster_result"
  result
}


#' Build the specificity-group network from a bluster_result
#'
#' @param bluster_result Output from \code{bluster_cluster()}.
#' @return An \code{igraph} object with node and edge attributes.
#' @export
bluster_network <- function(bluster_result) {

  stopifnot(inherits(bluster_result, "bluster_result"))

  g <- bluster_result$graph

  # Add cluster membership as vertex attribute
  mem <- bluster_result$membership
  cluster_map <- stats::setNames(mem$cluster_id, mem$cell_id)

  igraph::V(g)$cluster <- cluster_map[igraph::V(g)$name]

  # Add edge type attribute
  e_df <- bluster_result$edges
  igraph::E(g)$edge_type <- e_df$edge_type[seq_len(igraph::ecount(g))]

  g
}


# ---- internal scoring helpers ----------------------------------------

#' Collapse cells from the same clone to a representative
#' @keywords internal
.collapse_clones <- function(dt) {

  # Keep one representative per clone per chain type
  has_clone <- dt[!is.na(clone_id) & clone_id != ""]
  no_clone  <- dt[is.na(clone_id) | clone_id == ""]

  if (nrow(has_clone) > 0) {
    representatives <- has_clone[, .SD[1], by = .(clone_id, chain_type)]
    dt <- rbind(representatives, no_clone, fill = TRUE)
  }

  dt
}

#' Build combined edges from motifs and global results
#' @keywords internal
.build_combined_edges <- function(motif_results, global_edges, bcr_data,
                                  n_cores = 1L) {

  edge_list <- list()

  # Convert motif results to pairwise edges
  if (nrow(motif_results) > 0) {
    build_one <- function(i) {
      members <- motif_results$member_ids[[i]]
      if (length(members) < 2) return(NULL)

      # Create edges between all pairs sharing this motif
      pairs <- utils::combn(members, 2)
      data.table::data.table(
        seq_id_1  = pairs[1, ],
        seq_id_2  = pairs[2, ],
        weight    = 1.0,
        edge_type = "local_motif",
        motif     = motif_results$motif[i]
      )
    }

    motif_edges_list <- .bluster_lapply(seq_len(nrow(motif_results)),
                                        build_one, n_cores)
    motif_edges_list <- motif_edges_list[
      !vapply(motif_edges_list, is.null, logical(1))]
    edge_list <- c(edge_list, motif_edges_list)
  }

  # Add global similarity edges
  if (nrow(global_edges) > 0) {
    global_e <- data.table::data.table(
      seq_id_1  = global_edges$seq_id_1,
      seq_id_2  = global_edges$seq_id_2,
      weight    = 1.0 / (1.0 + global_edges$distance),
      edge_type = "global_similarity",
      motif     = NA_character_
    )
    edge_list[[length(edge_list) + 1]] <- global_e
  }

  if (length(edge_list) == 0) {
    return(data.table::data.table(
      seq_id_1 = character(0), seq_id_2 = character(0),
      weight = numeric(0), edge_type = character(0),
      motif = character(0)
    ))
  }

  edges <- data.table::rbindlist(edge_list, fill = TRUE)
  # Deduplicate (keep strongest edge between any pair)
  edges[, pair_key := paste(pmin(seq_id_1, seq_id_2),
                            pmax(seq_id_1, seq_id_2), sep = "::")]
  edges <- edges[, .SD[which.max(weight)], by = pair_key]
  edges[, pair_key := NULL]

  edges
}

#' Test V-gene enrichment in a cluster via Fisher's exact test
#' @keywords internal
.test_vgene_enrichment <- function(cluster_vgenes, all_vgenes) {

  cluster_vgenes <- cluster_vgenes[!is.na(cluster_vgenes)]
  all_vgenes <- all_vgenes[!is.na(all_vgenes)]

  if (length(cluster_vgenes) < 2 || length(all_vgenes) < 10) {
    return(list(pvalue = 1.0, gene = NA_character_))
  }

  # Find most common V gene in cluster
  tab <- sort(table(cluster_vgenes), decreasing = TRUE)
  top_gene <- names(tab)[1]

  # 2×2 contingency table
  a <- tab[1]                           # cluster, has gene
 b <- length(cluster_vgenes) - a       # cluster, other genes
  c_count <- sum(all_vgenes == top_gene) - a  # background, has gene
  d <- length(all_vgenes) - a - b - c_count   # background, other genes

  mat <- matrix(c(a, b, max(c_count, 0), max(d, 0)), nrow = 2)

  pval <- tryCatch(
    stats::fisher.test(mat, alternative = "greater")$p.value,
    error = function(e) 1.0
  )

  list(pvalue = pval, gene = top_gene)
}

#' Compute consensus CDR3 sequence for a cluster
#' @keywords internal
.consensus_cdr3 <- function(cdr3_vec) {

  cdr3_vec <- cdr3_vec[!is.na(cdr3_vec)]
  if (length(cdr3_vec) == 0) return(NA_character_)
  if (length(cdr3_vec) == 1) return(cdr3_vec)

  # Use most common CDR3 as representative
  tab <- sort(table(cdr3_vec), decreasing = TRUE)
  names(tab)[1]
}

#' Get the most frequent gene
#' @keywords internal
.top_gene <- function(genes) {
  genes <- genes[!is.na(genes)]
  if (length(genes) == 0) return(NA_character_)
  tab <- sort(table(genes), decreasing = TRUE)
  names(tab)[1]
}

#' Compute composite blusteR score
#' @keywords internal
.compute_bluster_score <- function(n_members, clonal_diversity,
                                  sample_diversity, vh_pvalue, len_sd) {

  # Higher is better
  # Rewards: many members, high clonal diversity, cross-sample convergence,
  #          significant V-gene enrichment, CDR3 length homogeneity
  size_score    <- log2(n_members + 1) / 10
  clone_score   <- clonal_diversity
  sample_score  <- sample_diversity
  vgene_score   <- pmin(-log10(vh_pvalue + 1e-300) / 50, 1.0)
  length_score  <- pmax(1 - len_sd / 5, 0)

  composite <- (0.20 * size_score +
                0.25 * clone_score +
                0.25 * sample_score +
                0.15 * vgene_score +
                0.15 * length_score)

  round(composite, 4)
}

#' Empty result constructor
#' @keywords internal
.empty_bluster_result <- function(bcr_data, motif_results) {
  result <- list(
    clusters   = data.table::data.table(),
    membership = data.table::data.table(),
    edges      = data.table::data.table(),
    graph      = igraph::make_empty_graph(),
    motifs     = motif_results,
    params     = list()
  )
  class(result) <- "bluster_result"
  result
}
