#!/usr/bin/env python3
"""Bar chart of the most enriched GO / pathway terms among the DE genes.

This is the phenotype-facing panel: each bar is a biological process that is
over-represented among the genes separating the tolerant and susceptible lines,
ranked by significance and coloured by fold enrichment. Hover shows the
overlapping genes.

Styled in the naraen.net "Phalaena Automata" mauve palette; with --png (and
kaleido) it also writes a static image.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go


PAL = {
    "bg":    "#1D171A",  # coffee-bean
    "grid":  "#392D34",  # shadow-grey
    "text":  "#D7C9D0",  # pale-slate
}
# Low-to-high fold-enrichment ramp across the mauve palette.
COLORSCALE = [[0.0, "#5B4752"], [0.5, "#967386"], [1.0, "#C2ADB8"]]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--enrichment", required=True, type=Path)
    p.add_argument("--topn", type=int, default=15)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--png", type=Path, default=None)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.enrichment, sep="\t")
    df = df.sort_values("padj", kind="mergesort").head(args.topn).iloc[::-1]
    if df.empty:
        raise SystemExit("enrichment table is empty")

    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    labels = [f"{n}" for n in df["name"]]

    fig = go.Figure(go.Bar(
        x=df["neg_log10_padj"],
        y=labels,
        orientation="h",
        marker=dict(
            color=df["enrichment"],
            colorscale=COLORSCALE,
            colorbar=dict(title="fold<br>enrich.", tickfont=dict(color=PAL["text"])),
            line=dict(width=0),
        ),
        customdata=np.stack(
            [df["term"], df["k_overlap"], df["set_size"], df["padj"],
             df["overlap_genes"]], axis=-1),
        hovertemplate=(
            "<b>%{y}</b> (%{customdata[0]})<br>"
            "overlap = %{customdata[1]}/%{customdata[2]} genes<br>"
            "padj = %{customdata[3]:.2e}<br>"
            "genes: %{customdata[4]}<extra></extra>"
        ),
    ))
    fig.update_layout(
        title=dict(
            text="<b>Enriched processes among DE genes</b><br>"
                 "<sup>tolerant vs. susceptible — hypergeometric ORA, BH-adjusted</sup>",
            x=0.02, xanchor="left", font=dict(color=PAL["text"]),
        ),
        xaxis_title="-log<sub>10</sub>(adjusted p-value)",
        template="plotly_dark",
        paper_bgcolor=PAL["bg"], plot_bgcolor=PAL["bg"],
        font=dict(color=PAL["text"], family="ui-monospace, monospace"),
        height=max(360, 34 * len(df) + 150),
        margin=dict(l=240, r=40, t=90, b=60),
    )
    fig.update_xaxes(gridcolor=PAL["grid"], zerolinecolor=PAL["grid"])
    fig.update_yaxes(automargin=True, tickfont=dict(size=11))

    fig.write_html(args.out, include_plotlyjs="cdn", full_html=True,
                   config={"displaylogo": False, "responsive": True})
    print(f"wrote {args.out} ({len(df)} terms)")
    if args.png is not None:
        fig.write_image(args.png, width=1100, height=max(360, 34 * len(df) + 150),
                        scale=2)
        print(f"wrote {args.png}")


if __name__ == "__main__":
    main()
