"""Minimal tests for the numerical study. Run from the project root:
    python -m pytest numerics/tests -q
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "numerics"))

from src.simulate import GAMMA0, simulate_hmm, stationary_distribution
from src.hmm_em import fit_hmm, forward_backward, validate_fit
from src.mixture_hmm_em import fit_mixture_hmm, validate_mixture_fit


def test_stationary_distribution():
    pi = stationary_distribution(GAMMA0)
    assert np.allclose(pi @ GAMMA0, pi, atol=1e-10)
    assert np.isclose(pi.sum(), 1.0)


def test_simulation_shapes_and_values():
    y, z, refined = simulate_hmm(300, seed=1)
    assert y.shape == z.shape == refined.shape == (300,)
    assert set(np.unique(z)) <= {0, 1}
    assert set(np.unique(refined)) <= {0, 1, 2}
    # refined state 2 iff aggregated state 2
    assert np.all((refined == 2) == (z == 1))
    assert np.all(np.isfinite(y))


def test_forward_backward_posteriors_sum_to_one():
    y, _, _ = simulate_hmm(200, seed=2)
    pi = np.array([0.5, 0.5])
    A = np.array([[0.9, 0.1], [0.1, 0.9]])
    loglik, gamma, xi_sum = forward_backward(
        y, pi, A, np.array([-1.0, 3.0]), np.array([1.0, 1.0]))
    assert np.isfinite(loglik)
    assert np.allclose(gamma.sum(axis=1), 1.0, atol=1e-10)
    assert np.isclose(xi_sum.sum(), len(y) - 1, atol=1e-6)


def test_em_fit_valid_and_monotone():
    y, _, _ = simulate_hmm(400, seed=3)
    for K in (2, 3):
        fit = fit_hmm(y, K, seed=10 + K, n_starts=4, max_iter=150)
        validate_fit(fit, len(y))  # transitions, sds>0, posteriors, no NaN,
        #                            EM monotonicity (within tolerance)
        assert fit.K == K
        assert np.all(np.diff(fit.means) >= 0), "states must be sorted by mean"


def test_bic_penalises_parameters():
    y, _, _ = simulate_hmm(400, seed=4)
    fit = fit_hmm(y, 2, seed=5, n_starts=3)
    assert fit.bic(400) > -2 * fit.loglik


def test_mixture_em_fit_valid_and_monotone():
    y, _, _ = simulate_hmm(260, seed=6)
    fit = fit_mixture_hmm(y, 2, M=2, seed=7, n_starts=2, max_iter=120)
    validate_mixture_fit(fit, len(y))
    assert fit.K == 2
    assert fit.M == 2
    assert np.all(np.diff(fit.state_means()) >= 0), "states must be sorted by mean"


def test_screened_multistart_paths_are_valid():
    y, _, _ = simulate_hmm(240, seed=8)
    fit = fit_hmm(y, 3, seed=9, n_starts=5, max_iter=60,
                  screen_iter=8, refine_top=2)
    validate_fit(fit, len(y))
    mix = fit_mixture_hmm(y, 2, M=2, seed=10, n_starts=5, max_iter=60,
                          screen_iter=8, refine_top=2)
    validate_mixture_fit(mix, len(y))


def test_output_files_exist():
    out = ROOT / "numerics" / "output"
    fig = ROOT / "article" / "figures"
    expected = [
        out / "review_gaussian_fit_summary.csv",
        out / "review_mixture_fit_summary.csv",
        out / "review_model_comparison.csv",
        out / "review_replications.csv",
        out / "review_summary_by_T.csv",
        out / "review_metadata.json",
        out / "crosstab_fit2.csv",
        out / "crosstab_fit3.csv",
        out / "residual_acf_fit2.csv",
        out / "residual_acf_fit3.csv",
        fig / "simulation_density.pdf",
        fig / "simulation_density.png",
        fig / "state_decoding_heatmap.pdf",
        fig / "bic_boxplot.pdf",
        fig / "bic_by_T.pdf",
        fig / "posterior_corrective_state.pdf",
        fig / "residual_acf.pdf",
        ROOT / "article" / "tables" / "simulation_table.tex",
        ROOT / "article" / "tables" / "replication_summary_table.tex",
        ROOT / "supplement" / "tables" / "review_summary_table.tex",
    ]
    missing = [str(p) for p in expected if not p.exists()]
    assert not missing, f"missing output files: {missing}"
