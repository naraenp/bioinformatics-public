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
