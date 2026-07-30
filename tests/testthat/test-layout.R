context("layout")

data("nba_clusts")

test_that("non-core layout works", {
    expect_is(clustree(nba_clusts, prefix = "K", use_core_edges = FALSE),
              c("gg", "ggplot"))
})

test_that("highlighting core works", {
    expect_is(clustree(nba_clusts, prefix = "K", highlight_core = TRUE),
              c("gg", "ggplot"))
})

test_that("alternative layouts work", {
    expect_is(clustree(nba_clusts, prefix = "K", layout = "tree"),
              c("gg", "ggplot"))
    expect_is(clustree(nba_clusts, prefix = "K", layout = "sugiyama"),
              c("gg", "ggplot"))
})

test_that("unordered clustree layout works", {
    expect_is(clustree(nba_clusts, prefix = "K", order_clusters = FALSE),
              c("gg", "ggplot"))
})

test_that("resolution labels can be turned off", {
    expect_is(clustree(nba_clusts, prefix = "K", show_res_labels = FALSE),
              c("gg", "ggplot"))
})

test_that("clustree layout rows are evenly spaced and centred", {
    layout <- clustree(nba_clusts, prefix = "K", return = "layout")

    rows <- split(layout$x, layout$y)

    # Each row is centred on zero
    centres <- vapply(rows, function(x) mean(range(x)), numeric(1))
    expect_equal(unname(centres), rep(0, length(centres)))

    # Clusters within a row are one unit apart
    gaps <- unique(unlist(lapply(rows, function(x) round(diff(sort(x)), 8))))
    expect_equal(gaps, 1)

    # There is one row per clustering, with the lowest resolution at the top
    expect_equal(length(rows), 5)
    expect_equal(sort(unique(layout$y)), -4:0)
})

test_that("clusters are ordered by their dominant parent", {
    # Cluster 1 of the first clustering contributes 0.4 of itself to cluster 1
    # of the second clustering and 0.6 to cluster 2, so cluster 2 is placed
    # first. Cluster 3 belongs to cluster 2 of the first clustering.
    clusterings <- cbind(
        K1 = c(1, 1, 1, 1, 1, 2, 2, 2),
        K2 = c(2, 2, 2, 1, 1, 3, 3, 3)
    )

    ordered <- order_clustree_clusters(clusterings)

    expect_equal(ordered[["K1"]], c(1, 2))
    expect_equal(ordered[["K2"]], c(2, 1, 3))
})

test_that("cluster sorting is numeric when it can be", {
    expect_equal(sort_clusters(c(10, 2, 1)), c(1, 2, 10))
    expect_equal(sort_clusters(c("10", "2", "1")), c("1", "2", "10"))
    expect_equal(sort_clusters(c("b", "c", "a")), c("a", "b", "c"))
})
