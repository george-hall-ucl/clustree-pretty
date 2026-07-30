#' Sort cluster labels
#'
#' Sort a vector of cluster labels numerically when every label can be
#' interpreted as a number and lexicographically otherwise.
#'
#' @param clusters vector of cluster labels
#'
#' @return sorted vector of cluster labels
#'
#' @keywords internal
sort_clusters <- function(clusters) {

    clusters <- unique(clusters)

    is_num <- suppressWarnings(all(!is.na(as.numeric(as.character(clusters)))))

    if (is_num) {
        clusters[order(as.numeric(as.character(clusters)))]
    } else {
        clusters[order(as.character(clusters))]
    }
}


#' Order clusters for the clustree layout
#'
#' Order the clusters at each resolution so that each cluster is placed next to
#' the cluster at the previous resolution that it is most closely related to.
#' Working down the tree, each cluster is assigned to the parent that
#' contributes the greatest proportion of its own samples to it. The clusters
#' assigned to a parent are then placed together, ordered by the proportion of
#' the parent they represent. This keeps related clusters near each other and
#' greatly reduces the number of crossing edges.
#'
#' @param clusterings numeric matrix containing clustering information, each
#' column contains clustering at a separate resolution
#'
#' @return list with one element per resolution containing the ordered cluster
#' labels for that resolution
#'
#' @keywords internal
order_clustree_clusters <- function(clusterings) {

    res_names <- colnames(clusterings)
    n_res <- length(res_names)

    clusters_by_res <- lapply(res_names, function(res) {
        sort_clusters(clusterings[, res])
    })

    ordered <- vector("list", n_res)
    # The first resolution has no parents to order against
    ordered[[1]] <- clusters_by_res[[1]]

    for (idx in seq_len(n_res - 1)) {

        parents <- ordered[[idx]]
        children <- clusters_by_res[[idx + 1]]

        parents_chr <- as.character(parents)
        children_chr <- as.character(children)

        from_clust <- factor(as.character(clusterings[, res_names[idx]]),
                             levels = parents_chr)
        to_clust <- factor(as.character(clusterings[, res_names[idx + 1]]),
                           levels = children_chr)

        # Proportion of each parent's samples that move to each child. Rows are
        # parents (in their already ordered order), columns are children.
        trans <- table(from_clust, to_clust)
        trans <- trans / rowSums(trans)

        # The parent contributing the greatest proportion of itself to each
        # child, which is the parent that child is placed under
        dominant <- parents_chr[apply(trans, 2, which.max)]

        ordered[[idx + 1]] <- unlist(lapply(seq_along(parents), function(p_idx) {
            is_child <- dominant == parents_chr[p_idx]

            if (!any(is_child)) {
                return(NULL)
            }

            props <- trans[parents_chr[p_idx], children_chr[is_child]]

            children[is_child][order(props, decreasing = TRUE)]
        }))
    }

    names(ordered) <- res_names

    return(ordered)
}


#' Get clustree layout positions
#'
#' Calculate node positions for the clustree layout. Each resolution is placed
#' on its own row, with the lowest resolution at the top. Within a row clusters
#' are evenly spaced and centred, in the order given by
#' [order_clustree_clusters()]. Because spacing is constant, rows containing
#' more clusters are wider than those containing fewer.
#'
#' @param graph [tidygraph::tbl_graph] object containing the tree graph
#' @param clusterings numeric matrix containing clustering information, each
#' column contains clustering at a separate resolution
#' @param prefix string indicating columns containing clustering information
#' @param order_clusters logical, whether to order clusters within each
#' resolution using [order_clustree_clusters()]. If `FALSE` clusters are placed
#' in sorted order.
#'
#' @return data.frame with `x` and `y` columns giving the position of each node,
#' in the order the nodes appear in `graph`
#'
#' @keywords internal
clustree_layout_positions <- function(graph, clusterings, prefix,
                                      order_clusters = TRUE) {

    res_names <- colnames(clusterings)

    if (order_clusters) {
        ordered <- order_clustree_clusters(clusterings)
    } else {
        ordered <- lapply(res_names, function(res) {
            sort_clusters(clusterings[, res])
        })
    }

    positions <- lapply(seq_along(res_names), function(idx) {
        clusters <- ordered[[idx]]
        n_clusters <- length(clusters)
        res_clean <- as.numeric(gsub(prefix, "", res_names[idx]))

        data.frame(
            node = paste0(prefix, res_clean, "C", clusters),
            # Evenly spaced and centred on zero
            x = seq_len(n_clusters) - (n_clusters + 1) / 2,
            # Lowest resolution at the top
            y = -(idx - 1),
            stringsAsFactors = FALSE
        )
    })

    positions <- do.call("rbind", positions)

    node_names <- igraph::vertex_attr(graph, "node")
    positions <- positions[match(node_names, positions$node), ]

    return(positions)
}
