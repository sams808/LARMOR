# HMQC & correlation (2D)

For heteronuclear-correlation data. Open the `2rr` on the contour map.

## Overlay a 1D on a projection
Press **HMQC**, then **F2 1D…** / **F1 1D…**: click the spectrum in the left
**Explorer** (or Browse to one elsewhere). It is overlaid on that projection and
highlighted in the axis colour (F2 orange, F1 purple).

## Scale
The overlay is peak-matched automatically; press **fit scale** for a robust
least-squares scaling over the visible range, or type a manual multiplier.

## Un-correlated features
**uncorrelated F2 →** / **uncorrelated F1 →** sends `1D − scale·projection` to a
new workspace: the projection is the *correlated* signal, so the difference is
what does **not** show a cross-peak. Fit it like any 1D.

## Navigate
Extracting a trace opens a **new workspace**; the map stays in its own. Switch
via the **Workspaces** dock or **Back to 2D map** (Ctrl+2).

*(A generalized multi-experiment correlation — 1D + MQMAS + HMQC + REDOR — is on
the roadmap.)*
