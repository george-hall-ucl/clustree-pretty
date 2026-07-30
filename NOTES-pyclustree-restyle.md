# Restyling clustree after pyclustree — working notes

Handoff notes for continuing this work. Not part of the package (see
`.Rbuildignore`). Delete before opening a PR.

## Where things stand

Branch `claude/clustree-visualization-styling-7idu0x`, commit `e64d936`, pushed.
No PR opened.

**The next task is to redo this as a smaller, tighter PR.** The user's words:
"try not to add too much code. with each addition, consider whether it really
needs to be added. I don't want this code to be rejected as AI slop". So treat
`e64d936` as a *spike*: the algorithm in it is verified and worth keeping, the
surface area around it is not. Read "What to cut" below before writing anything.

## What the user actually asked for

> I like the clustree approach in this package but I think they are ugly. I much
> prefer the look of the clustrees in pyclustree, for example
> <https://pyclustree.complextissue.com/latest/_images/example_7_0.png>

Then, after a first attempt that only changed colours:

> these don't look anything like the ones from pyclustree ... keep rendering
> plots and seeing how the position of the nodes can be made more like
> pyclustree (**this is the main thing I want**)

Node positions are the deliverable. Everything else is negotiable.

## Getting the reference image

The docs site is blocked by the agent proxy (403 on CONNECT), so
`curl`/WebFetch on the `_images/` URL fails. The image is embedded as base64 in
the docs notebook, which *is* reachable:

```bash
curl -sSL -o example.ipynb \
  https://raw.githubusercontent.com/complextissue/pyclustree/main/docs/source/example.ipynb
python3 -c "
import json,base64
nb=json.load(open('example.ipynb'))
for i,c in enumerate(nb['cells']):
    for o in c.get('outputs',[]):
        if 'image/png' in o.get('data',{}):
            open(f'ref_cell{i}.png','wb').write(base64.b64decode(o['data']['image/png']))"
```

`ref_cell7.png` is exactly the image the user linked. Look at it before
changing anything — this whole task is visual, and reading it directly is worth
more than any description of it.

Source to read: `pyclustree/_clustree.py` and `pyclustree/_utils.py` on
raw.githubusercontent.com.

## The two things that create the pyclustree look

clustree had neither. Everything else (colours, edges, fonts) is cosmetic
detail on top.

1. **`nx.multipartite_layout(align="horizontal")`** — within a row, nodes are
   evenly spaced at *constant* spacing and centred on zero. Rows with more
   clusters are therefore wider, which is what produces the flask silhouette.
   clustree used igraph Reingold-Tilford, where x comes from subtree structure —
   that is the lopsided look the user objected to.

2. **`order_unique_clusters`** (`_utils.py`) — decides the left-to-right order
   within each row. This is what makes most edges short and vertical.

Precise semantics of (2), which are easy to get subtly wrong:

- Transition matrix is `pd.crosstab(parent, child, normalize="index")`, i.e.
  **row-normalised**: `T[p, c] = (cells going p→c) / (size of p)`. Note this is
  *not* clustree's `in_prop`, which is column-normalised
  (`count / size of child`). Using `in_prop` here gives different, worse orders.
- Each child is assigned to the parent maximising `T[p, c]` (`idxmax(axis=0)`).
- Parents are visited in their own already-ordered order; each parent's children
  are emitted sorted by `T[parent, child]` **descending**.
- Row 1 is plain sorted order (numeric sort when labels parse as numbers).

Positions: `x = seq_len(n) - (n + 1) / 2` (centred, unit spacing),
`y = -(level - 1)` (lowest resolution at top). That's all `multipartite_layout`
amounts to here; no need to replicate its rescaling, ggplot handles extents.

## How to verify the ordering (do this, don't eyeball it)

This was the highest-value step and it is cheap to repeat. Run pyclustree's real
Python function against the R implementation on identical data and diff.

```bash
pip install -q pandas          # not present by default
```

Export a clusterings matrix from R to CSV, run `order_unique_clusters` verbatim
(copy it out of `_utils.py`, it has no pyclustree dependencies — only pandas),
write one line per level, and `diff` against the R output. On a 14-resolution
tree this produced **byte-identical output**, which is what makes it safe to
claim the layout matches rather than approximates.

Scratch copies of both scripts were in the session scratchpad
(`py_order.py`, `export_clusts.R`); they are ~20 lines each, quicker to rewrite
than to recover.

Keep this check as a throwaway — it should *not* end up in the PR (it would add
a Python dependency to an R package).

## Environment setup

R is not installed in a fresh container. This works and takes a few minutes:

```bash
apt-get update -qq
apt-get install -y -qq r-base-core r-cran-ggplot2 r-cran-igraph r-cran-dplyr \
  r-cran-rlang r-cran-checkmate r-cran-viridis r-cran-ggrepel r-cran-tidygraph \
  r-cran-ggraph r-cran-devtools
```

All deps are available as Debian packages — do **not** try `install.packages()`,
it is much slower. Seurat and SingleCellExperiment are not available; the ~16
tests that need them skip cleanly, which is expected and not a failure.

Rendering loop:

```r
suppressMessages(devtools::load_all("."))
data(nba_clusts)
ggplot2::ggsave("out.png", clustree(nba_clusts, prefix = "K"),
                width = 8, height = 7, dpi = 150)
```

then Read the PNG. `nba_clusts` only has 5 resolutions — **too small to judge
the silhouette**. Build a bigger tree (14 k-means clusterings at k = 5..12 over
12 gaussian blobs) and render at `width=10, height=15` to compare against
`ref_cell7.png`. Without that you cannot tell whether the layout is right.

## What is in `e64d936`, and what to cut

Roughly 217 new lines of R plus 150 changed in `clustree.R`. Honest split:

**Keep — this is the actual ask.**

- `R/layout.R` (150 lines): `sort_clusters`, `order_clustree_clusters`,
  `clustree_layout_positions`. This is the verified algorithm. Could lose maybe
  20 lines of roxygen but the logic is about as tight as it gets.
- In `clustree.matrix`: `layout = c("clustree", "tree", "sugiyama")` default plus
  the `create_layout(graph, "manual", ...)` branch. ~6 lines.
- Panel padding (`expansion(mult = 0.05, add = ...)` on x and y). This is a
  genuine bug fix, not polish: node size is in points but the panel only
  expanded in data units, so nodes at the edges of the tree were clipped. It is
  visible on a 2-resolution tree. Keep, but it can be shorter than what is there
  now.

**Cut — dead code.**

- `clustree_edge_colour` in `R/theme.R`: defined, documented, exported to
  `man/`, and **never referenced** — the default is hardcoded `"grey35"` in the
  `clustree.matrix` signature. Verified with grep. This is exactly the kind of
  thing that reads as AI slop. Delete it and `man/clustree_edge_colour.Rd`.

**Cut or defer — separate features riding along with the layout change.**

- `show_res_labels` / `res_label_size` and the row-label block (~25 lines plus
  the `left_pad` heuristic). It is genuinely what allows dropping the redundant
  colour legend, and it matches the reference, but it is a second feature. If
  the goal is a reviewable PR, this belongs in its own.
- `R/clustree_overlay.R` changes (palette consistency, 15 lines). Scope creep —
  the complaint was about `clustree()`, not `clustree_overlay()`.
- `edge_colour` argument + `col2rgb` validation, `highlight_core` rework,
  closed arrowheads, bold white node text. Each is 2–10 lines and defensible,
  but together they are why the diff got big. Consider shipping the layout
  alone first.
- `add_clustree_colour_scale` + viridis fallback pulls two new `NAMESPACE`
  imports. If the categorical palette stays, this can probably collapse into a
  few lines inline rather than a helper with its own Rd page.

**Note on scale:** 8 new `man/*.Rd` files are generated for internal helpers.
That is correct per the repo's existing convention (`build_tree_graph.Rd` etc.
exist), but it inflates the diff. Fewer helpers means fewer Rd files.

## Judgement calls made, still open

- **Kept the node-size legend**, which the reference does not have. It encodes
  information nothing else shows. It is also the last thing making the plot
  asymmetric. `+ guides(size = "none")` removes it. User has not weighed in.
- **Changed defaults**, so existing users' plots change. Intended, given the
  ask, but it is a breaking visual change for a CRAN package. `NEWS.md` has an
  entry under "development version"; `DESCRIPTION` version deliberately left
  alone.
- `clustree_palette_n(nrow(nodes))` in the overlay — `nrow(nodes)` is an upper
  bound on discrete levels. Replaced a magic `100`, which would have errored
  above 100 levels.

## Housekeeping that must not regress

- `devtools::document()` bumps `RoxygenNote` in `DESCRIPTION` and rewrites
  `man/clustree-package.Rd` with author/URL boilerplate. Both are unrelated
  churn — revert them every time:
  `git checkout -- man/clustree-package.Rd DESCRIPTION`
- The vignette had two chunks using edge aesthetics that no longer exist
  (`scale_edge_color_continuous`, `guides(edge_colour = FALSE)`). If edge
  colouring is reworked again, `vignettes/clustree.Rmd` needs updating or it
  documents a broken API. Check by purling and evaluating every chunk, not by
  reading it.
- Full suite: `testthat::test_dir("tests/testthat")` after `load_all()`.
  49 pass, ~16 skip (Seurat/SCE absent).

## Suggested next steps

1. Branch fresh from `main`. Do not build on `e64d936`.
2. Port `R/layout.R` across, plus the `layout` argument and the padding fix.
   Aim for the smallest diff that changes node positions.
3. Re-run the Python parity check to confirm the port is still faithful.
4. Render `nba_clusts` and a 14-resolution tree; compare against
   `ref_cell7.png`.
5. Keep the layout tests (`test-layout.R`); they are cheap and cover the
   ordering rule with a hand-computed case.
6. Decide with the user whether row labels / colours / edge styling go in this
   PR or a follow-up, rather than deciding for them.
