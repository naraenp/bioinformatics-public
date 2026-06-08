#!/usr/bin/env python3
"""Render an interactive Plotly volcano plot from DE results.

Genes are colored by significance (up / down / ns) under the supplied FDR
and |log2FC| thresholds. A curated set of canonical AML marker genes is
annotated by name so the plot reads as a quick sanity check.

With --png and kaleido installed, also writes a static PNG (used for the
portfolio thumbnail). Styled in the naraen.net "Phalaena Automata" mauve
palette on a dark ground.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go


AML_SIGNATURE = [
    "FLT3", "KIT", "MYC", "BCL2", "MEIS1", "HOXA9", "WT1", "CD34", "MPO",
    "CEBPA", "RUNX1", "DNMT3A", "TET2", "IDH1", "IDH2", "NPM1", "ASXL1",
    "EZH2", "GATA2", "SPI1",
]


# Phalaena Automata palette (see naraen.net).
PAL = {
    "bg":        "#1D171A",  # coffee-bean
    "grid":      "#392D34",  # shadow-grey
    "ns":        "#5B4752",  # mauve-shadow
    "down":      "#967386",  # dusty-mauve (accent)
    "up":        "#C2ADB8",  # lilac-ash
    "marker":    "#EAE1E5",  # alabaster (signature outline + labels)
    "text":      "#D7C9D0",  # pale-slate
    "thresh":    "#5B4752",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--de", required=True, type=Path)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--lfc", type=float, default=1.0)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--png", type=Path, default=None,
                   help="optional static PNG export (needs kaleido)")
    return p.parse_args()


def classify(df: pd.DataFrame, fdr: float, lfc: float) -> pd.Series:
    sig = df["padj"] < fdr
    up = sig & (df["log2FC"] >= lfc)
    down = sig & (df["log2FC"] <= -lfc)
    label = np.where(up, "Up in AML", np.where(down, "Down in AML", "Not significant"))
    return pd.Series(label, index=df.index, name="status")


def build_figure(df: pd.DataFrame, fdr: float, lfc: float) -> go.Figure:
    df = df.copy()
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["status"] = classify(df, fdr, lfc)

    color_map = {
        "Up in AML":       PAL["up"],
        "Down in AML":     PAL["down"],
        "Not significant": PAL["ns"],
    }

    fig = go.Figure()
    for status in ["Not significant", "Down in AML", "Up in AML"]:
        sub = df[df["status"] == status]
        if sub.empty:
            continue
        fig.add_trace(go.Scattergl(
            x=sub["log2FC"],
            y=sub["neg_log10_padj"],
            mode="markers",
            name=f"{status} (n={len(sub):,})",
            marker=dict(
                color=color_map[status],
                size=6,
                opacity=0.55 if status == "Not significant" else 0.85,
                line=dict(width=0),
            ),
            customdata=np.stack([sub["gene"], sub["padj"], sub["pvalue"]], axis=-1),
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>"
                "log2FC = %{x:.2f}<br>"
                "-log10(padj) = %{y:.2f}<br>"
                "padj = %{customdata[1]:.2e}<br>"
                "p = %{customdata[2]:.2e}<extra></extra>"
            ),
        ))

    # Threshold guide lines.
    fig.add_hline(y=-np.log10(fdr), line=dict(color=PAL["thresh"], dash="dash", width=1))
    fig.add_vline(x=lfc,  line=dict(color=PAL["thresh"], dash="dash", width=1))
    fig.add_vline(x=-lfc, line=dict(color=PAL["thresh"], dash="dash", width=1))

    # Annotate the curated AML marker set.
    sig_df = df[df["gene"].isin(AML_SIGNATURE)].copy()
    if not sig_df.empty:
        fig.add_trace(go.Scatter(
            x=sig_df["log2FC"],
            y=sig_df["neg_log10_padj"],
            mode="markers+text",
            text=sig_df["gene"],
            textposition="top center",
            textfont=dict(size=11, color=PAL["marker"]),
            marker=dict(
                color="rgba(0,0,0,0)",
                size=12,
                line=dict(color=PAL["marker"], width=1.2),
            ),
            name="AML markers",
            hoverinfo="skip",
            showlegend=True,
        ))

    fig.update_layout(
        title=dict(
            text="<b>AML vs. healthy bulk RNA-seq</b><br>"
                 "<sup>Welch t-test on log2-CPM, BH FDR; "
                 f"thresholds: |log2FC| ≥ {lfc}, FDR &lt; {fdr}</sup>",
            x=0.02, xanchor="left", font=dict(color=PAL["text"]),
        ),
        xaxis_title="log<sub>2</sub> fold change (AML / healthy)",
        yaxis_title="-log<sub>10</sub>(adjusted p-value)",
        template="plotly_dark",
        paper_bgcolor=PAL["bg"],
        plot_bgcolor=PAL["bg"],
        font=dict(color=PAL["text"], family="ui-monospace, monospace"),
        height=720,
        margin=dict(l=70, r=30, t=90, b=60),
        legend=dict(
            orientation="h",
            yanchor="bottom", y=1.02,
            xanchor="right",  x=1.0,
            bgcolor="rgba(29,23,26,0.6)",
        ),
        hoverlabel=dict(bgcolor=PAL["bg"], font_size=12),
    )
    fig.update_xaxes(zeroline=True, zerolinecolor=PAL["grid"], gridcolor=PAL["grid"])
    fig.update_yaxes(zerolinecolor=PAL["grid"], gridcolor=PAL["grid"])
    return fig


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.de, sep="\t")
    fig = build_figure(df, args.fdr, args.lfc)
    fig.write_html(
        args.out,
        include_plotlyjs="cdn",
        full_html=True,
        config={"displaylogo": False, "responsive": True},
    )
    print(f"wrote {args.out} ({(df['padj'] < args.fdr).sum()} genes at FDR<{args.fdr})")

    if args.png is not None:
        fig.write_image(args.png, width=1200, height=720, scale=2)
        print(f"wrote {args.png}")


if __name__ == "__main__":
    main()
