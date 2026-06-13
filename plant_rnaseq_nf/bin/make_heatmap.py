#!/usr/bin/env python3
"""Clustered heatmap of the top differentially expressed genes.

Takes the per-gene normalized expression (log2 normalized counts) and the DE
table, picks the strongest DE genes, row-z-scores them, orders samples by
genotype, and renders an interactive Plotly heatmap. Genes are ordered by
hierarchical clustering when SciPy is available, otherwise by log2 fold change.

Styled in the naraen.net "Aequorea" abyssal-marine palette so it sits natively
in the portfolio. With --png (and kaleido) it also writes a static image.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go


# Aequorea palette (see naraen.net / make_volcano.py): abyssal marine + glow.
PAL = {
    "bg":     "#0A141A",  # abyss
    "grid":   "#163239",  # fathom
    "low":    "#4E7C86",  # muted teal   (down in tolerant)
    "mid":    "#0A141A",  # abyss        (z = 0)
    "high":   "#2F9FB0",  # teal         (up in tolerant)
    "text":   "#BBD7DC",  # haze
}
# Diverging, zero-centred colourscale: soft mist on the negative arm sinking to
# the abyssal centre, climbing to a bioluminescent glow on the positive arm.
COLORSCALE = [
    [0.0, "#8CB6BE"], [0.25, PAL["low"]], [0.5, PAL["mid"]],
    [0.75, PAL["high"]], [1.0, "#4BDDE6"],
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--norm", required=True, type=Path, help="norm_counts.tsv")
    p.add_argument("--de", required=True, type=Path, help="de_results.tsv")
    p.add_argument("--meta", required=True, type=Path, help="metadata.tsv")
    p.add_argument("--topn", type=int, default=40)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--lfc", type=float, default=1.0)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--png", type=Path, default=None)
    return p.parse_args()


def cluster_order(mat: np.ndarray) -> np.ndarray:
    """Leaf order from hierarchical clustering; identity order on failure."""
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage
        if mat.shape[0] < 3:
            return np.arange(mat.shape[0])
        return leaves_list(linkage(mat, method="average", metric="correlation"))
    except Exception:
        return np.arange(mat.shape[0])


def main() -> None:
    args = parse_args()
    norm = pd.read_csv(args.norm, sep="\t", index_col=0)
    de = pd.read_csv(args.de, sep="\t").dropna(subset=["padj"])
    meta = pd.read_csv(args.meta, sep="\t")

    sig = de[(de["padj"] < args.fdr) & (de["log2FC"].abs() >= args.lfc)]
    chosen = (sig if not sig.empty else de).head(args.topn)
    genes = [g for g in chosen["gene"] if g in norm.index]
    if not genes:
        raise SystemExit("no DE genes found in the normalized matrix")

    # Order samples by genotype for a readable block structure.
    meta = meta.sort_values(["genotype", "sample_id"])
    samples = [s for s in meta["sample_id"].astype(str) if s in norm.columns]
    sub = norm.loc[genes, samples]

    # Row z-score (expression relative to the gene's own mean).
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, 1), axis=0)
    row_order = cluster_order(z.to_numpy())
    z = z.iloc[row_order]

    geno = meta.set_index("sample_id")["genotype"].astype(str)
    col_labels = [f"{s}<br>({geno.get(s, '?')})" for s in samples]

    fig = go.Figure(go.Heatmap(
        z=z.to_numpy(),
        x=col_labels,
        y=z.index.tolist(),
        colorscale=COLORSCALE,
        zmid=0,
        colorbar=dict(title="row z", tickfont=dict(color=PAL["text"])),
        hovertemplate="gene %{y}<br>%{x}<br>z = %{z:.2f}<extra></extra>",
    ))
    fig.update_layout(
        title=dict(
            text=f"<b>Top {len(genes)} differentially expressed genes</b><br>"
                 "<sup>tolerant vs. susceptible — log2 normalized counts, "
                 "row z-scored</sup>",
            x=0.02, xanchor="left", font=dict(color=PAL["text"]),
        ),
        template="plotly_dark",
        paper_bgcolor=PAL["bg"], plot_bgcolor=PAL["bg"],
        font=dict(color=PAL["text"], family="ui-monospace, monospace"),
        height=max(480, 18 * len(genes) + 160),
        margin=dict(l=120, r=40, t=90, b=90),
    )
    fig.update_xaxes(tickfont=dict(size=9))
    fig.update_yaxes(tickfont=dict(size=10))

    fig.write_html(args.out, include_plotlyjs="cdn", full_html=True,
                   config={"displaylogo": False, "responsive": True})
    print(f"wrote {args.out} ({len(genes)} genes x {len(samples)} samples)")
    if args.png is not None:
        fig.write_image(args.png, width=1100,
                        height=max(480, 18 * len(genes) + 160), scale=2)
        print(f"wrote {args.png}")


if __name__ == "__main__":
    main()
