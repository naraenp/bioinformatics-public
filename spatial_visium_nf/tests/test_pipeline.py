"""Fast, mostly data-free unit tests for the spatial_visium_nf helpers.

These cover the pure logic — NNLS deconvolution, Moran's I + the kNN graph,
per-cell-type signature means, and the synthetic-demo generators — without
invoking the full scanpy/squidpy stack on real data; that path is exercised by
`run_local.sh --demo`. Run from the repo root:  pytest spatial_visium_nf/tests/
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


# ---- deconvolve_nnls.nnls_deconvolve ---------------------------------------

def test_nnls_recovers_known_mixture():
    pytest.importorskip("scipy")
    pytest.importorskip("pandas")
    import numpy as np
    dec = load("deconvolve_nnls")

    rng = np.random.default_rng(0)
    G, K, n = 8, 3, 25
    S = rng.uniform(0.1, 1.0, size=(G, K))           # genes x types
    P_true = rng.dirichlet(np.ones(K), size=n)        # spots x types, rows sum 1
    X = P_true @ S.T                                  # exact linear mixtures
    P = dec.nnls_deconvolve(X, S)
    assert P.shape == (n, K)
    assert np.allclose(P.sum(axis=1), 1.0)
    assert np.allclose(P, P_true, atol=1e-6)


def test_gene_scale_factors_mean_and_floor():
    import numpy as np
    dec = load("deconvolve_nnls")
    S = np.array([[2.0, 4.0],      # mean 3
                  [0.0, 0.0],      # mean 0 -> floored
                  [1.0, 1.0]])     # mean 1
    g = dec.gene_scale_factors(S, floor=1e-6)
    assert np.allclose(g, [3.0, 1e-6, 1.0])
    assert (g > 0).all()                              # safe to divide by


def test_gene_scaled_nnls_still_recovers_known_mixture():
    pytest.importorskip("scipy")
    import numpy as np
    dec = load("deconvolve_nnls")
    rng = np.random.default_rng(7)
    G, K, n = 10, 3, 20
    S = rng.uniform(0.1, 5.0, size=(G, K))            # wide per-gene scale spread
    P_true = rng.dirichlet(np.ones(K), size=n)
    X = P_true @ S.T
    g = dec.gene_scale_factors(S)
    P = dec.nnls_deconvolve(X / g[None, :], S / g[:, None])
    # exact linear mixtures are recovered regardless of the per-gene weighting
    assert np.allclose(P, P_true, atol=1e-6)


def test_nnls_proportions_nonnegative_and_normalised():
    pytest.importorskip("scipy")
    pytest.importorskip("pandas")
    import numpy as np
    dec = load("deconvolve_nnls")
    rng = np.random.default_rng(1)
    X = rng.uniform(0, 1, size=(10, 6))               # noisy, not in span
    S = rng.uniform(0, 1, size=(6, 4))
    P = dec.nnls_deconvolve(X, S)
    assert (P >= 0).all()
    assert np.allclose(P.sum(axis=1), 1.0)


# ---- svg_moran: kNN weights + Moran's I ------------------------------------

def test_knn_weights_symmetric_zero_diagonal():
    pytest.importorskip("scipy")
    import numpy as np
    svg = load("svg_moran")
    coords = np.array([[x, y] for x in range(5) for y in range(5)], float)
    W = svg.knn_weights(coords, k=4)
    assert np.allclose(W, W.T)                         # symmetric
    assert np.allclose(np.diag(W), 0.0)               # no self-loops
    assert (W.sum(axis=1) >= 4).all()                 # at least k neighbours


def test_morans_i_high_on_gradient_low_on_noise():
    pytest.importorskip("scipy")
    import numpy as np
    svg = load("svg_moran")
    coords = np.array([[x, y] for x in range(8) for y in range(8)], float)
    W = svg.knn_weights(coords, k=4)
    rng = np.random.default_rng(2)
    gradient = coords[:, 0]                            # smooth spatial trend
    noise = rng.normal(size=coords.shape[0])
    V = np.column_stack([gradient, noise])
    I = svg.morans_i(V, W)
    assert I[0] > 0.5                                  # gradient: strong autocorr
    assert abs(I[1]) < 0.3                             # noise: near zero


def test_permutation_pvals_bounded_and_separates_signal():
    pytest.importorskip("scipy")
    import numpy as np
    svg = load("svg_moran")
    coords = np.array([[x, y] for x in range(8) for y in range(8)], float)
    W = svg.knn_weights(coords, k=4)
    rng = np.random.default_rng(2)
    V = np.column_stack([coords[:, 0], rng.normal(size=coords.shape[0])])
    I = svg.morans_i(V, W)
    n_perm = 200
    p = svg.permutation_pvals(V, W, I, n_perm=n_perm, seed=0)
    # One-sided permutation p-values live in (0, 1]; the smallest attainable is
    # 1/(n_perm+1) because the observed statistic is counted in the numerator.
    assert np.all(p > 0) and np.all(p <= 1.0)
    assert p.min() >= 1.0 / (n_perm + 1) - 1e-12
    # The smooth gradient is autocorrelated, the noise gene is not.
    assert p[0] < p[1]
    assert p[0] == pytest.approx(1.0 / (n_perm + 1))  # gradient hits the floor
    # Seeded, so the test is deterministic across runs.
    assert np.array_equal(p, svg.permutation_pvals(V, W, I, n_perm=n_perm, seed=0))


# ---- build_signature.mean_by_label -----------------------------------------

def test_mean_by_label():
    pytest.importorskip("pandas")
    import numpy as np
    bs = load("build_signature")
    X = np.array([[1.0, 0.0], [3.0, 0.0], [0.0, 4.0], [0.0, 6.0]])
    labels = ["A", "A", "B", "B"]
    M, types = bs.mean_by_label(X, labels)
    assert types == ["A", "B"]
    assert M.shape == (2, 2)                           # genes x types
    assert np.allclose(M[:, 0], [2.0, 0.0])           # mean of the A cells
    assert np.allclose(M[:, 1], [0.0, 5.0])           # mean of the B cells


# ---- make_demo_data pure-numpy cores ---------------------------------------

def test_signatures_and_proportions_are_distributions():
    import numpy as np
    dd = load("make_demo_data")
    rng = np.random.default_rng(3)
    S, marker_idx = dd.build_signatures(rng, n_genes=60, n_types=4,
                                        markers_per_type=5, marker_fold=8.0)
    assert S.shape == (60, 4)
    assert np.allclose(S.sum(axis=0), 1.0)            # each signature a profile
    # a type's marker block is elevated in its own column vs the mean column
    for k, idx in marker_idx.items():
        assert S[idx, k].mean() > S[idx].mean()

    coords = [(r, c) for r in range(6) for c in range(6)]
    P = dd.planted_proportions(rng, coords, n_types=4, sigma=2.0, concentration=12.0)
    assert P.shape == (36, 4)
    assert (P >= 0).all()
    assert np.allclose(P.sum(axis=1), 1.0)

    rate = dd.mix_spot_expression(P, S)
    assert np.allclose(rate.sum(axis=1), 1.0)


def test_demo_data_generation(tmp_path):
    pytest.importorskip("h5py")
    pytest.importorskip("scipy")
    import sys
    dd = load("make_demo_data")
    argv = ["make_demo_data", "--outdir", str(tmp_path), "--grid", "6",
            "--cells-per-type", "20", "--genes", "80", "--seed", "1"]
    old = sys.argv
    sys.argv = argv
    try:
        dd.main()
    finally:
        sys.argv = old
    assert (tmp_path / "spatial" / "filtered_feature_bc_matrix.h5").exists()
    assert (tmp_path / "spatial" / "spatial" / "tissue_positions_list.csv").exists()
    assert (tmp_path / "reference" / "matrix.mtx.gz").exists()
    assert (tmp_path / "reference" / "metadata.csv").exists()
    prop = (tmp_path / "truth" / "proportions.csv").read_text().splitlines()
    assert prop[0].startswith("barcode,")
    assert len(prop) == 1 + 36                        # 6x6 spots + header
