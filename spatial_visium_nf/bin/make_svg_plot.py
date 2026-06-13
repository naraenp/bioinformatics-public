#!/usr/bin/env python3
"""MAKE_SVG_PLOT — interactive spatial expression of the top spatially variable genes.

Takes the Moran's I ranking and the spatial expression matrix and renders a
dropdown of the top spatially variable genes; selecting a gene shows its
log-normalised expression across the spots in array space. Marine "Aequorea"
palette; self-contained, iframe-ready HTML.

Input:  svg.tsv + spatial_norm.h5ad (layers['lognorm'], obsm['spatial'])
Output: spatial_svg.html
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

PAL = {"bg": "#0A141A", "grid": "#163239", "text": "#BBD7DC"}  # Aequorea: abyss / fathom / haze
EXPR_SCALE = [[0.0, "#0A141A"], [0.35, "#294E57"], [0.7, "#4BDDE6"], [1.0, "#DBEBEE"]]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--svg", required=True, type=Path, help="svg.tsv")
    p.add_argument("--spatial", required=True, type=Path, help="spatial_norm.h5ad")
    p.add_argument("--topn", type=int, default=6)
    p.add_argument("--out", required=True, type=Path, help="spatial_svg.html")
    p.add_argument("--png", type=Path, default=None)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    import scanpy as sc

    svg = pd.read_csv(args.svg, sep="\t")
    adata = sc.read_h5ad(args.spatial)
    coords = np.asarray(adata.obsm["spatial"], float)
    x, y = coords[:, 0], coords[:, 1]
    X = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)

    genes = [g for g in svg["gene"].head(args.topn) if g in set(adata.var_names)]
    if not genes:
        raise SystemExit("no top SVG genes found in the spatial matrix")
    gidx = {g: i for i, g in enumerate(adata.var_names)}
    imap = svg.set_index("gene")["morans_I"]

    fig = go.Figure()
    for n, g in enumerate(genes):
        expr = X[:, gidx[g]]
        fig.add_trace(go.Scatter(
            x=x, y=y, mode="markers", name=g, visible=(n == 0),
            marker=dict(size=9, color=expr, colorscale=EXPR_SCALE,
                        colorbar=dict(title="log1p", tickfont=dict(color=PAL["text"])),
                        line=dict(width=0.3, color=PAL["grid"])),
            hovertemplate=f"{g} = %{{marker.color:.2f}}<br>spot %{{text}}<extra></extra>",
            text=adata.obs_names,
        ))

    buttons = []
    for n, g in enumerate(genes):
        vis = [i == n for i in range(len(genes))]
        buttons.append(dict(label=f"{g}  (I={imap.get(g, float('nan')):.2f})",
                            method="update",
                            args=[{"visible": vis},
                                  {"title.text": f"<b>Spatial expression — {g}</b>"
                                   f"<br><sup>Moran's I = {imap.get(g, float('nan')):.3f}"
                                   " — top spatially variable genes</sup>"}]))

    g0 = genes[0]
    fig.update_layout(
        title=dict(text=f"<b>Spatial expression — {g0}</b><br>"
                        f"<sup>Moran's I = {imap.get(g0, float('nan')):.3f}"
                        " — top spatially variable genes</sup>",
                   x=0.02, xanchor="left", font=dict(color=PAL["text"])),
        template="plotly_dark",
        paper_bgcolor=PAL["bg"], plot_bgcolor=PAL["bg"],
        font=dict(color=PAL["text"], family="ui-monospace, monospace"),
        updatemenus=[dict(buttons=buttons, x=1.0, xanchor="right", y=1.16,
                          yanchor="top", bgcolor=PAL["grid"],
                          font=dict(color=PAL["text"]))],
        showlegend=False, height=620, margin=dict(l=50, r=40, t=110, b=50),
    )
    fig.update_xaxes(title="array_col", showgrid=False, zeroline=False,
                     scaleanchor="y", scaleratio=1)
    fig.update_yaxes(title="array_row", autorange="reversed", showgrid=False,
                     zeroline=False)

    fig.write_html(args.out, include_plotlyjs="cdn", full_html=True,
                   config={"displaylogo": False, "responsive": True})
    print(f"wrote {args.out} (top {len(genes)} SVGs)")
    if args.png is not None:
        try:                                  # static export is best-effort (needs kaleido)
            fig.write_image(args.png, width=1000, height=620, scale=2)
            print(f"wrote {args.png}")
        except Exception as e:                # noqa: BLE001 — don't fail the stage on it
            print(f"WARNING: static PNG export skipped ({type(e).__name__}); "
                  f"the interactive {args.out} is the canonical output")


if __name__ == "__main__":
    main()
