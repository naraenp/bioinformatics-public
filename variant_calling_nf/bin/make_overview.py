#!/usr/bin/env python3
"""Cohort variant-overview chart: per-sample variant counts by genic region.

Reads the annotated.tsv from annotate_variants.py and renders an interactive
grouped bar chart — for each sample, how many PASS variants it carries
(dosage > 0) split into exonic / intronic / intergenic — with the cohort's
SNV transition/transversion ratio in the subtitle (a standard germline-callset
QC read-out). Styled in the naraen.net "Aequorea" abyssal-marine palette so it
sits natively in the portfolio.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

# Aequorea palette (see naraen.net / sibling pipelines): abyssal marine + glow.
PAL = {"bg": "#0A141A", "grid": "#163239", "text": "#BBD7DC"}
REGION_COLORS = {
    "exonic": "#4BDDE6",       # bioluminescent glow — the consequential ones
    "intronic": "#2F9FB0",     # teal
    "intergenic": "#4E7C86",   # muted teal
}
REGIONS = ["exonic", "intronic", "intergenic"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--annotated", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--png", type=Path, default=None)
    return p.parse_args()


def tstv_ratio(df: pd.DataFrame) -> float | None:
    snv = df[df["kind"] == "snv"]
    ts = int((snv["tstv"] == "ts").sum())
    tv = int((snv["tstv"] == "tv").sum())
    return ts / tv if tv else None


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.annotated, sep="\t")
    df = df[df["filter"].isin([".", "PASS"])]

    meta_cols = ["chrom", "pos", "ref", "alt", "kind", "tstv", "region", "filter"]
    samples = [c for c in df.columns if c not in meta_cols]

    # Count, per sample, variants present (dosage > 0) in each region.
    counts = {r: [] for r in REGIONS}
    for s in samples:
        dose = pd.to_numeric(df[s], errors="coerce").fillna(0)
        present = df[dose > 0]
        vc = present["region"].value_counts()
        for r in REGIONS:
            counts[r].append(int(vc.get(r, 0)))

    fig = go.Figure()
    for r in REGIONS:
        fig.add_bar(name=r, x=samples, y=counts[r],
                    marker_color=REGION_COLORS[r])

    ratio = tstv_ratio(df)
    sub = (f"per-sample PASS variants by genic region — cohort Ts/Tv = {ratio:.2f}"
           if ratio is not None else "per-sample PASS variants by genic region")
    fig.update_layout(
        barmode="group",
        title=dict(text=f"<b>Cohort variant overview</b><br><sup>{sub}</sup>",
                   x=0.02, xanchor="left", font=dict(color=PAL["text"])),
        template="plotly_dark",
        paper_bgcolor=PAL["bg"], plot_bgcolor=PAL["bg"],
        font=dict(color=PAL["text"], family="ui-monospace, monospace"),
        legend=dict(title="region", bgcolor="rgba(0,0,0,0)"),
        height=480, margin=dict(l=70, r=40, t=90, b=60),
    )
    fig.update_xaxes(title="sample", gridcolor=PAL["grid"])
    fig.update_yaxes(title="variant count", gridcolor=PAL["grid"])

    fig.write_html(args.out, include_plotlyjs="cdn", full_html=True,
                   config={"displaylogo": False, "responsive": True})
    print(f"wrote {args.out} ({len(samples)} samples)")
    if args.png is not None:
        fig.write_image(args.png, width=1100, height=480, scale=2)
        print(f"wrote {args.png}")


if __name__ == "__main__":
    main()
