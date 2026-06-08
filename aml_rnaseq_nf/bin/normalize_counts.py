#!/usr/bin/env python3
"""CPM normalize a raw count matrix and return log2(CPM + 1).

Uses per-sample library sizes; not a full TMM, but adequate for a teaching
pipeline and stable for the volcano downstream.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    counts = pd.read_csv(args.counts, sep="\t", index_col=0)
    lib_sizes = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.divide(lib_sizes, axis=1) * 1e6
    lcpm = np.log2(cpm.fillna(0) + 1.0)
    lcpm.to_csv(args.out, sep="\t")
    print(f"wrote {args.out} (shape={lcpm.shape}, median lib size={int(lib_sizes.median())})")


if __name__ == "__main__":
    main()
