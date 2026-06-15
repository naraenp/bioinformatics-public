#!/usr/bin/env python3
"""Genotype heatmap of the most discriminating variants across the cohort.

Reads annotated.tsv (per-sample alt-allele dosage 0/1/2) and renders an
interactive Plotly heatmap of the top-variance variants — samples on x, variant
sites on y, colour = dosage. Samples are ordered by group (from the
samplesheet) so any group structure shows as a clean block, the variant-calling
analog of the genotype-block heatmap in plant_rnaseq_nf. Styled in the
naraen.net "Aequorea" abyssal-marine palette.

With a single-sample callset (e.g. the real GIAB run) it still renders — one
column of the carried variants.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

PAL = {"bg": "#0A141A", "grid": "#163239", "text": "#BBD7DC"}
# 0 = abyss, 1 = teal, 2 = bioluminescent glow.
COLORSCALE = [[0.0, "#0A141A"], [0.5, "#2F9FB0"], [1.0, "#4BDDE6"]]
META_COLS = ["chrom", "pos", "ref", "alt", "kind", "tstv", "region", "filter"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--annotated", required=True, type=Path)
    p.add_argument("--meta", type=Path, default=None,
                   help="samplesheet.csv (sample_id, group) for column ordering")
    p.add_argument("--topn", type=int, default=40)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--png", type=Path, default=None)
    return p.parse_args()


def order_samples(samples: list[str], meta: Path | None) -> tuple[list[str], dict]:
    if meta is None or not meta.exists():
        return samples, {}
    m = pd.read_csv(meta)
    grp = dict(zip(m["sample_id"].astype(str), m["group"].astype(str)))
    ordered = sorted(samples, key=lambda s: (grp.get(s, "?"), s))
    return ordered, grp


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.annotated, sep="\t")
    df = df[df["filter"].isin([".", "PASS"])].reset_index(drop=True)
    samples = [c for c in df.columns if c not in META_COLS]
    if not samples:
        raise SystemExit("no sample columns in annotated.tsv")

    dose = df[samples].apply(pd.to_numeric, errors="coerce")
    # Rank variants by cross-sample variance: the discriminating sites float up.
    var = dose.var(axis=1).fillna(0)
    keep = var.sort_values(ascending=False).head(args.topn).index
    dose = dose.loc[keep]
    labels = [f"{df.loc[i,'chrom']}:{df.loc[i,'pos']} {df.loc[i,'ref']}>{df.loc[i,'alt']}"
              for i in keep]

    samples, grp = order_samples(samples, args.meta)
    mat = dose[samples].to_numpy(dtype=float)
    col_labels = [f"{s}<br>({grp[s]})" if s in grp else s for s in samples]

    fig = go.Figure(go.Heatmap(
        z=mat, x=col_labels, y=labels,
        colorscale=COLORSCALE, zmin=0, zmax=2,
        colorbar=dict(title="alt dosage", tickvals=[0, 1, 2],
                      tickfont=dict(color=PAL["text"])),
        hovertemplate="%{y}<br>%{x}<br>dosage = %{z}<extra></extra>",
    ))
    fig.update_layout(
        title=dict(
            text=f"<b>Top {len(labels)} discriminating variants</b>"
                 "<br><sup>alt-allele dosage per sample; sites ranked by "
                 "cross-sample variance</sup>",
            x=0.02, xanchor="left", font=dict(color=PAL["text"])),
        template="plotly_dark",
        paper_bgcolor=PAL["bg"], plot_bgcolor=PAL["bg"],
        font=dict(color=PAL["text"], family="ui-monospace, monospace"),
        height=max(420, 16 * len(labels) + 160),
        margin=dict(l=170, r=40, t=90, b=80),
    )
    fig.update_xaxes(tickfont=dict(size=9))
    fig.update_yaxes(tickfont=dict(size=9))

    fig.write_html(args.out, include_plotlyjs="cdn", full_html=True,
                   config={"displaylogo": False, "responsive": True})
    print(f"wrote {args.out} ({len(labels)} variants x {len(samples)} samples)")
    if args.png is not None:
        fig.write_image(args.png, width=1000,
                        height=max(420, 16 * len(labels) + 160), scale=2)
        print(f"wrote {args.png}")


if __name__ == "__main__":
    main()
