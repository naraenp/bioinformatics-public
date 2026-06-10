#!/usr/bin/env python3
"""MAKE_CELLTYPE_PLOT — interactive spatial map of the deconvolution.

Renders the spots in array space with a dropdown that switches between the
dominant cell type per spot (categorical) and each cell type's per-spot
proportion (continuous mauve map) — the spatial-omics analogue of the plant
pipeline's heatmap, and a preview of the Part B viewer. Styled in the naraen.net
"Phalaena Automata" mauve palette and emitted as a self-contained, iframe-ready
HTML.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

# Phalaena Automata palette (see naraen.net / plant_rnaseq_nf make_heatmap.py).
PAL = {"bg": "#1D171A", "grid": "#392D34", "text": "#D7C9D0",
       "low": "#1D171A", "high": "#C2ADB8"}
PROP_SCALE = [[0.0, "#1D171A"], [0.35, "#5B4752"], [0.7, "#967386"], [1.0, "#EAE1E5"]]
# Categorical colours for dominant-type view (mauve-family, distinguishable).
CAT = ["#C2ADB8", "#967386", "#7D9B8E", "#B8946F", "#8F8AB0", "#A0707F", "#6E8FA8"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--proportions", required=True, type=Path, help="proportions.tsv")
    p.add_argument("--out", required=True, type=Path, help="spatial_celltypes.html")
    p.add_argument("--png", type=Path, default=None)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.proportions, sep="\t", index_col=0)
    meta_cols = {"dominant_type", "array_row", "array_col"}
    types = [c for c in df.columns if c not in meta_cols]
    if "array_col" in df and "array_row" in df:
        x, y = df["array_col"], df["array_row"]
    else:
        raise SystemExit("proportions.tsv lacks array_row/array_col coordinates")

    fig = go.Figure()
    # Trace 0..K-1: dominant-type, one scatter per type (categorical legend).
    for j, t in enumerate(types):
        m = df["dominant_type"] == t
        fig.add_trace(go.Scatter(
            x=x[m], y=y[m], mode="markers", name=t, visible=True,
            marker=dict(size=9, color=CAT[j % len(CAT)],
                        line=dict(width=0.5, color=PAL["grid"])),
            hovertemplate=f"{t}<br>spot %{{text}}<extra></extra>",
            text=df.index[m],
        ))
    # Trace K..2K-1: per-type continuous proportion maps (hidden initially).
    for t in types:
        fig.add_trace(go.Scatter(
            x=x, y=y, mode="markers", name=f"{t} prop", visible=False,
            marker=dict(size=9, color=df[t], colorscale=PROP_SCALE, cmin=0, cmax=1,
                        colorbar=dict(title="prop", tickfont=dict(color=PAL["text"])),
                        line=dict(width=0.3, color=PAL["grid"])),
            hovertemplate=f"{t} = %{{marker.color:.2f}}<br>spot %{{text}}<extra></extra>",
            text=df.index,
        ))

    n = len(types)
    def vis(dominant: bool, k: int | None = None):
        v = [False] * (2 * n)
        if dominant:
            for i in range(n):
                v[i] = True
        else:
            v[n + k] = True
        return v

    buttons = [dict(label="Dominant type", method="update",
                    args=[{"visible": vis(True)},
                          {"showlegend": True}])]
    for k, t in enumerate(types):
        buttons.append(dict(label=f"{t} proportion", method="update",
                            args=[{"visible": vis(False, k)},
                                  {"showlegend": False}]))

    fig.update_layout(
        title=dict(text="<b>Per-spot cell-type deconvolution (NNLS)</b><br>"
                        "<sup>dominant type, or pick a cell type's proportion map</sup>",
                   x=0.02, xanchor="left", font=dict(color=PAL["text"])),
        template="plotly_dark",
        paper_bgcolor=PAL["bg"], plot_bgcolor=PAL["bg"],
        font=dict(color=PAL["text"], family="ui-monospace, monospace"),
        updatemenus=[dict(buttons=buttons, x=1.0, xanchor="right", y=1.16,
                          yanchor="top", bgcolor=PAL["grid"],
                          font=dict(color=PAL["text"]))],
        legend=dict(title="dominant", bgcolor="rgba(0,0,0,0)"),
        height=620, margin=dict(l=50, r=40, t=110, b=50),
    )
    # Visium array orientation: row increases downward.
    fig.update_xaxes(title="array_col", showgrid=False, zeroline=False,
                     scaleanchor="y", scaleratio=1)
    fig.update_yaxes(title="array_row", autorange="reversed", showgrid=False,
                     zeroline=False)

    fig.write_html(args.out, include_plotlyjs="cdn", full_html=True,
                   config={"displaylogo": False, "responsive": True})
    print(f"wrote {args.out} ({df.shape[0]} spots, {n} cell types)")
    if args.png is not None:
        try:                                  # static export is best-effort (needs kaleido)
            fig.write_image(args.png, width=1000, height=620, scale=2)
            print(f"wrote {args.png}")
        except Exception as e:                # noqa: BLE001 — don't fail the stage on it
            print(f"WARNING: static PNG export skipped ({type(e).__name__}); "
                  f"the interactive {args.out} is the canonical output")


if __name__ == "__main__":
    main()
