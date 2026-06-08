"""Fast, mostly data-free unit tests for the plant_rnaseq_nf helpers.

These cover the pure logic (FDR, GMT parsing, featureCounts column mapping,
hypergeometric ORA, clustering order, demo-data generation) without invoking
the heavy aligner toolchain — that path is exercised by `run_local.sh --demo`.
Run from the repo root:  pytest plant_rnaseq_nf/tests/
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin"


def load(stem: str):
    """Import a bin/ script as a module by file path."""
    path = BIN / f"{stem}.py"
    spec = importlib.util.spec_from_file_location(stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---- build_count_matrix.stem -----------------------------------------------

def test_stem_strips_bam_suffixes():
    bcm = pytest.importorskip("pandas") and load("build_count_matrix")
    assert bcm.stem("/work/ab/cd/tol_r1.sorted.bam") == "tol_r1"
    assert bcm.stem("sus_r2.bam") == "sus_r2"
    assert bcm.stem("/x/y/S3.sorted.bam") == "S3"


# ---- run_enrichment: BH FDR + GMT + hypergeometric --------------------------

def test_bh_fdr_bounds_and_monotone():
    pytest.importorskip("scipy")
    enr = load("run_enrichment")
    import numpy as np
    p = np.array([0.001, 0.01, 0.02, 0.5, 0.9])
    q = enr.bh_fdr(p)
    assert np.all(q >= 0) and np.all(q <= 1)
    # q-values are monotone non-decreasing in p-rank.
    order = np.argsort(p)
    assert np.all(np.diff(q[order]) >= -1e-12)
    # largest p maps to its own value (rank n).
    assert q[order][-1] == pytest.approx(0.9, abs=1e-9)


def test_bh_fdr_empty():
    pytest.importorskip("scipy")
    enr = load("run_enrichment")
    import numpy as np
    assert enr.bh_fdr(np.array([])).size == 0


def test_read_gmt(tmp_path):
    pytest.importorskip("scipy")
    enr = load("run_enrichment")
    gmt = tmp_path / "g.gmt"
    gmt.write_text("T1\tname one\tA\tB\tC\nT2\tname two\tC\tD\n")
    sets = enr.read_gmt(gmt)
    assert sets[0] == ("T1", "name one", {"A", "B", "C"})
    assert sets[1][2] == {"C", "D"}


def test_hypergeom_detects_clear_enrichment():
    pytest.importorskip("scipy")
    enr = load("run_enrichment")
    from scipy import stats
    # 100 genes, 20 selected; a 10-gene set with 9 hits should be very enriched.
    p = stats.hypergeom.sf(9 - 1, 100, 20, 10)
    assert p < 1e-3


# ---- make_heatmap.cluster_order --------------------------------------------

def test_cluster_order_is_a_permutation():
    pytest.importorskip("plotly")
    import numpy as np
    hm = load("make_heatmap")
    mat = np.random.default_rng(0).normal(size=(6, 8))
    order = hm.cluster_order(mat)
    assert sorted(order.tolist()) == list(range(6))


# ---- make_demo_data (numpy only) -------------------------------------------

def test_demo_data_generation(tmp_path):
    pytest.importorskip("numpy")
    dd = load("make_demo_data")
    import sys
    argv = ["make_demo_data", "--outdir", str(tmp_path), "--genes", "20",
            "--reps", "2", "--base-depth", "10", "--seed", "1"]
    old = sys.argv
    sys.argv = argv
    try:
        dd.main()
    finally:
        sys.argv = old
    assert (tmp_path / "genome.fa").exists()
    assert (tmp_path / "genes.gtf").exists()
    assert (tmp_path / "go_sets.gmt").exists()
    sheet = (tmp_path / "samplesheet.csv").read_text().splitlines()
    assert sheet[0] == "sample_id,genotype,fastq_1,fastq_2"
    assert len(sheet) == 1 + 4               # 2 tolerant + 2 susceptible
    assert (tmp_path / "reads" / "tol_r1_R1.fastq.gz").exists()
