#' clustree colour palette
#'
#' A categorical colour palette (matplotlib's "tab20", as used by packages
#' such as `pyclustree`) used as the default discrete node colour scale for
#' clustering tree plots. A qualitative palette of muted, evenly spaced hues
#' is easier to visually separate than `ggplot2`'s default (highly
#' saturated) hue palette.
#'
#' @keywords internal
clustree_palette <- c(
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a",
    "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
    "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d",
    "#17becf", "#9edae5"
)

#' Default flat edge colour for clustree plots
#'
#' @keywords internal
clustree_edge_colour <- "grey35"

#' Get `n` colours from the clustree palette
#'
#' Recycles [clustree_palette] if there are more levels than colours.
#'
#' @param n number of colours required
#'
#' @keywords internal
clustree_palette_n <- function(n) {
    rep_len(clustree_palette, n)
}

#' Add a clustree colour/fill scale
#'
#' Adds a discrete scale using the clustree palette or, for numeric data
#' (e.g. gene expression), a continuous viridis scale, matching the
#' aesthetic used for that data.
#'
#' @param gg ggplot object to add the scale to
#' @param values vector of values being mapped to the aesthetic
#' @param aesthetic either "colour" or "fill"
#'
#' @keywords internal
#'
#' @importFrom ggplot2 scale_colour_manual scale_fill_manual
#' @importFrom viridis scale_color_viridis scale_fill_viridis
add_clustree_colour_scale <- function(gg, values, aesthetic = c("colour", "fill")) {

    aesthetic <- match.arg(aesthetic)

    if (is.numeric(values)) {
        if (aesthetic == "colour") {
            gg <- gg + viridis::scale_color_viridis()
        } else {
            gg <- gg + viridis::scale_fill_viridis()
        }
    } else {
        pal <- clustree_palette_n(length(unique(values)))
        if (aesthetic == "colour") {
            gg <- gg + scale_colour_manual(values = pal, drop = FALSE)
        } else {
            gg <- gg + scale_fill_manual(values = pal, drop = FALSE)
        }
    }

    gg
}
