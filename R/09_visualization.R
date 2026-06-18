#' Plot blusteR cluster network
#'
#' Visualises the specificity-group network.  Nodes are cells/clones,
#' coloured by cluster.  Edges are coloured by type (local motif vs
#' global similarity).
#'
#' @param bluster_result A \code{bluster_result} object.
#' @param layout igraph layout function (default
#'   \code{igraph::layout_with_fr}).
#' @param node_size Node size (default 3).
#' @param label_clusters Label cluster IDs (default TRUE for ≤ 30
#'   clusters).
#' @param highlight_annotated Highlight clusters with known
#'   specificity (default TRUE).
#'
#' @return A ggplot2 object.
#' @export
plot_cluster_network <- function(bluster_result,
                                 layout = NULL,
                                 node_size = 3,
                                 label_clusters = NULL,
                                 highlight_annotated = TRUE) {

  stopifnot(inherits(bluster_result, "bluster_result"))

  g <- bluster_result$graph
  if (igraph::vcount(g) == 0) {
    message("[blusteR] Empty graph; nothing to plot.")
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::ggtitle("No clusters to display"))
  }

  # Get layout coordinates
  if (is.null(layout)) {
    layout_fn <- igraph::layout_with_fr
  } else {
    layout_fn <- layout
  }
  coords <- layout_fn(g)
  colnames(coords) <- c("x", "y")

  # Build node data
  mem <- bluster_result$membership
  cluster_map <- stats::setNames(mem$cluster_id, mem$cell_id)

  node_df <- data.frame(
    name      = igraph::V(g)$name,
    x         = coords[, 1],
    y         = coords[, 2],
    cluster   = cluster_map[igraph::V(g)$name],
    stringsAsFactors = FALSE
  )
  node_df$cluster[is.na(node_df$cluster)] <- "unassigned"

  # Highlight annotated clusters
  if (highlight_annotated && !is.null(bluster_result$annotations) &&
      nrow(bluster_result$annotations) > 0) {
    ann_clusters <- unique(bluster_result$annotations$cluster_id)
    node_df$annotated <- node_df$cluster %in% ann_clusters
  } else {
    node_df$annotated <- FALSE
  }

  # Build edge data
  el <- igraph::as_edgelist(g)
  edge_df <- data.frame(
    x    = coords[match(el[,1], igraph::V(g)$name), 1],
    y    = coords[match(el[,1], igraph::V(g)$name), 2],
    xend = coords[match(el[,2], igraph::V(g)$name), 1],
    yend = coords[match(el[,2], igraph::V(g)$name), 2],
    stringsAsFactors = FALSE
  )

  if (!is.null(bluster_result$edges) && nrow(bluster_result$edges) > 0 &&
      "edge_type" %in% names(bluster_result$edges)) {
    edge_df$type <- bluster_result$edges$edge_type[seq_len(nrow(edge_df))]
  } else {
    edge_df$type <- "unknown"
  }

  # Auto-decide labelling
  n_cl <- length(unique(node_df$cluster[node_df$cluster != "unassigned"]))
  if (is.null(label_clusters)) label_clusters <- n_cl <= 30

  # Plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = type),
      alpha = 0.3, linewidth = 0.3
    ) +
    ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(x = x, y = y, fill = cluster),
      shape = 21, size = node_size, stroke = 0.3,
      colour = ifelse(node_df$annotated, "red", "grey30")
    ) +
    ggplot2::scale_colour_manual(
      values = c(local_motif = "#2166AC", global_similarity = "#B2182B",
                 unknown = "grey50"),
      name = "Edge type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text  = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = "blusteR Specificity-Group Network",
                  fill = "Cluster")

  # Add cluster labels at centroids
  if (label_clusters) {
    centroids <- do.call(rbind, lapply(
      unique(node_df$cluster[node_df$cluster != "unassigned"]),
      function(cl) {
        sub <- node_df[node_df$cluster == cl, ]
        data.frame(x = mean(sub$x), y = mean(sub$y),
                   label = cl, stringsAsFactors = FALSE)
      }
    ))
    if (!is.null(centroids) && nrow(centroids) > 0) {
      p <- p + ggplot2::geom_text(
        data = centroids,
        ggplot2::aes(x = x, y = y, label = label),
        size = 2.5, fontface = "bold"
      )
    }
  }

  p
}


#' Plot V-gene usage across clusters
#'
#' @param bluster_result A \code{bluster_result} object.
#' @param top_n Show the top N V genes (default 15).
#'
#' @return A ggplot2 object.
#' @export
plot_vgene_usage <- function(bluster_result, top_n = 15L) {

  cl <- bluster_result$clusters
  if (nrow(cl) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  # Extract V gene per cluster
  vg_data <- cl[, .(cluster_id, top_vh)]
  vg_data <- vg_data[!is.na(top_vh)]

  # Count
  vg_counts <- vg_data[, .N, by = top_vh]
  data.table::setorder(vg_counts, -N)
  vg_counts <- head(vg_counts, top_n)

  ggplot2::ggplot(vg_counts, ggplot2::aes(
    x = stats::reorder(top_vh, N), y = N
  )) +
    ggplot2::geom_col(fill = "#2166AC", alpha = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Top V-Gene Usage Across blusteR Clusters",
      x = "V Gene", y = "Number of Clusters"
    )
}


#' Plot motif enrichment results
#'
#' Displays a motif logo-like representation by showing the most
#' enriched motifs and their fold-enrichment over background.
#'
#' @param bluster_result A \code{bluster_result} object.
#' @param top_n Number of top motifs (default 20).
#'
#' @return A ggplot2 object.
#' @export
plot_motif_logo <- function(bluster_result, top_n = 20L) {

  motifs <- bluster_result$motifs
  if (is.null(motifs) || nrow(motifs) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::ggtitle("No enriched motifs"))
  }

  data.table::setorder(motifs, pvalue_adj)
  show <- head(motifs, top_n)

  show[, neg_log_p := -log10(pvalue_adj + 1e-300)]

  ggplot2::ggplot(show, ggplot2::aes(
    x = stats::reorder(motif, neg_log_p), y = neg_log_p
  )) +
    ggplot2::geom_col(ggplot2::aes(fill = chain), alpha = 0.85) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(heavy = "#2166AC", light = "#D6604D")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Top Enriched CDR3 Motifs",
      x = "Motif",
      y = expression(-log[10](p[adj])),
      fill = "Chain"
    )
}


#' Plot cluster summary dashboard
#'
#' Four-panel overview: cluster size distribution, blusteR score
#' distribution, clonal diversity vs size, and sample diversity.
#'
#' @param bluster_result A \code{bluster_result} object.
#'
#' @return A ggplot2 object (uses patchwork-style manual layout via
#'   gridExtra if available, otherwise returns size histogram).
#' @export
plot_cluster_summary <- function(bluster_result) {

  cl <- bluster_result$clusters
  if (nrow(cl) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  # Panel 1: Cluster size distribution
  p1 <- ggplot2::ggplot(cl, ggplot2::aes(x = n_members)) +
    ggplot2::geom_histogram(bins = 30, fill = "#2166AC", alpha = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Cluster Size Distribution", x = "Members", y = "Count")

  # Panel 2: blusteR score
  p2 <- ggplot2::ggplot(cl, ggplot2::aes(x = bluster_score)) +
    ggplot2::geom_histogram(bins = 30, fill = "#B2182B", alpha = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "blusteR Score Distribution", x = "Score", y = "Count")

  # Panel 3: Clonal diversity vs size
  p3 <- ggplot2::ggplot(cl, ggplot2::aes(
    x = n_members, y = clonal_diversity
  )) +
    ggplot2::geom_point(alpha = 0.6, colour = "#2166AC") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Clonal Diversity vs Size",
                  x = "Cluster Size", y = "Clonal Diversity")

  # Panel 4: Sample diversity
  p4 <- ggplot2::ggplot(cl, ggplot2::aes(x = n_samples)) +
    ggplot2::geom_bar(fill = "#4DAF4A", alpha = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Sample Representation", x = "Samples per Cluster",
                  y = "Count")

  # Try to combine with gridExtra
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  } else {
    message("Install 'gridExtra' for combined panel plot. Returning size histogram.")
    return(p1)
  }
}
