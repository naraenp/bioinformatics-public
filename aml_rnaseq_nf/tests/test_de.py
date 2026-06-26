"""Fast, data-free unit tests for the DE + normalization math.

These exercise the pure functions in bin/ on small synthetic inputs, so the
suite runs in well under a second and needs no recount3 download — only the
same conda env as the pipeline. Run from the repo root:

    pytest aml_rnaseq_nf/tests/
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_MISSING = [m for m in ("numpy", "pandas", "scipy")
            if importlib.util.find_spec(m) is None]
if _MISSING:
    pytest.skip(
        f"tests need the aml_rnaseq env (missing: {', '.join(_MISSING)})",
        allow_module_level=True,
    )

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

BIN = Path(__file__).resolve().parents[1] / "bin"


def _load(name: str):
    """Import a bin/ script by path as a module."""
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_bh_fdr_is_monotone_and_bounded() -> None:
    run_de = _load("run_de")
    pvals = np.array([0.001, 0.008, 0.02, 0.2, 0.5, 0.9])
    q = run_de.bh_fdr(pvals)
    assert np.all((q >= 0) & (q <= 1))
    # BH q-values are non-decreasing in the p-value rank.
    assert np.all(np.diff(q[np.argsort(pvals)]) >= -1e-12)
    # q >= p for every gene.
    assert np.all(q + 1e-12 >= pvals)


def test_bh_fdr_all_significant() -> None:
    run_de = _load("run_de")
    q = run_de.bh_fdr(np.full(50, 1e-6))
    assert np.all(q < 0.05)


def test_cpm_filter_drops_low_expression_genes() -> None:
    load_counts = _load("load_counts")
    counts = pd.DataFrame(
        {"s1": [1000, 0, 5000], "s2": [1200, 0, 4800], "s3": [900, 0, 5200], "s4": [1100, 0, 5100]},
        index=["EXPRESSED_A", "SILENT", "EXPRESSED_B"],
    )
    kept = load_counts.cpm_filter(counts, min_cpm=1.0, min_frac=0.25)
    assert "SILENT" not in kept.index
    assert {"EXPRESSED_A", "EXPRESSED_B"}.issubset(kept.index)


def test_volcano_classify_respects_fdr_and_lfc_thresholds() -> None:
    # make_volcano imports plotly at module top; the lean CI env omits it, so
    # skip there (the classify logic itself is pure pandas/numpy).
    pytest.importorskip("plotly")
    make_volcano = _load("make_volcano")
    df = pd.DataFrame({
        "gene":   ["UP", "DOWN", "NS_PVAL", "NS_LFC", "EDGE_UP", "EDGE_DOWN"],
        "padj":   [0.001, 0.001, 0.50, 0.001, 0.001, 0.001],
        "log2FC": [2.0, -2.0, 3.0, 0.5, 1.0, -1.0],
    })
    status = make_volcano.classify(df, fdr=0.05, lfc=1.0)
    by_gene = dict(zip(df["gene"], status))
    assert by_gene["UP"] == "Up in AML"
    assert by_gene["DOWN"] == "Down in AML"
    # Significant fold change but padj above the FDR cutoff -> not significant.
    assert by_gene["NS_PVAL"] == "Not significant"
    # Significant padj but |log2FC| below the threshold -> not significant.
    assert by_gene["NS_LFC"] == "Not significant"
    # Boundary: |log2FC| == lfc is inclusive (the code uses >= / <=).
    assert by_gene["EDGE_UP"] == "Up in AML"
    assert by_gene["EDGE_DOWN"] == "Down in AML"
