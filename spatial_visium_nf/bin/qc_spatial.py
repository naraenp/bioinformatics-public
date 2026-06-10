#!/usr/bin/env python3
"""QC_SPATIAL — per-spot quality control and filtering.

Computes total counts, genes-per-spot and mitochondrial fraction for each
Visium spot, drops spots failing the thresholds, and writes a small JSON QC
summary alongside the filtered spatial_qc.h5ad. Mitochondrial genes are detected
by name prefix (MT- for human; absent in the demo, which yields pct_mt = 0).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import scanpy as sc


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="inp", required=True, type=Path, help="spatial.h5ad")
    p.add_argument("--out", required=True, type=Path, help="spatial_qc.h5ad")
    p.add_argument("--qc-json", type=Path, default=None, help="QC summary JSON")
    p.add_argument("--min-counts", type=int, default=200)
    p.add_argument("--min-genes", type=int, default=50)
    p.add_argument("--max-pct-mt", type=float, default=30.0)
    p.add_argument("--mt-prefix", default="MT-")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    adata = sc.read_h5ad(args.inp)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(args.mt_prefix.upper())
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)

    n0 = adata.n_obs
    keep = (
        (adata.obs["total_counts"] >= args.min_counts)
        & (adata.obs["n_genes_by_counts"] >= args.min_genes)
        & (adata.obs["pct_counts_mt"] <= args.max_pct_mt)
    )
    adata = adata[keep].copy()

    summary = {
        "spots_in": int(n0),
        "spots_kept": int(adata.n_obs),
        "spots_dropped": int(n0 - adata.n_obs),
        "median_counts": float(np.median(adata.obs["total_counts"])),
        "median_genes": float(np.median(adata.obs["n_genes_by_counts"])),
        "max_pct_mt_kept": float(adata.obs["pct_counts_mt"].max()),
        "thresholds": {"min_counts": args.min_counts, "min_genes": args.min_genes,
                       "max_pct_mt": args.max_pct_mt},
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    if args.qc_json is not None:
        args.qc_json.parent.mkdir(parents=True, exist_ok=True)
        args.qc_json.write_text(json.dumps(summary, indent=2))
    print(f"wrote {args.out}: kept {adata.n_obs}/{n0} spots "
          f"(dropped {n0 - adata.n_obs})")


if __name__ == "__main__":
    main()
