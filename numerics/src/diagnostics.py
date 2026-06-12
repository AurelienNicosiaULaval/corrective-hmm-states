"""Diagnostics: fit summaries, model comparison, true-vs-fitted crosstabs."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .hmm_em import HMMFit
from .mixture_hmm_em import MixtureHMMFit
from .simulate import REFINED_LABELS


def fit_summary_rows(fit: HMMFit, T: int) -> list[dict]:
    """One row per fitted state: mean, sd, posterior occupation, logLik, BIC, ICL."""
    occ = fit.gamma.mean(axis=0)
    return [
        {
            "K": fit.K,
            "state": j + 1,
            "mean": fit.means[j],
            "sd": fit.sds[j],
            "posterior_occupation": occ[j],
            "logLik": fit.loglik,
            "BIC": fit.bic(T),
            "ICL": fit.icl(T),
        }
        for j in range(fit.K)
    ]


def crosstab_true_vs_fitted(refined: np.ndarray, fit: HMMFit,
                            use: str = "viterbi",
                            proportions: bool = False) -> pd.DataFrame:
    """Crosstab of true refined states {1a, 1b, 2} vs fitted state labels.

    use : "viterbi" (decoded path) or "posterior" (argmax of gamma).
    """
    fitted = fit.viterbi if use == "viterbi" else fit.gamma.argmax(axis=1)
    tab = pd.crosstab(
        pd.Series(REFINED_LABELS[refined], name="true_refined"),
        pd.Series(fitted + 1, name=f"fitted_state_K{fit.K}"),
    )
    if proportions:
        tab = tab.div(tab.sum(axis=1), axis=0)
    return tab


def model_comparison(fits: dict[int, HMMFit], T: int) -> dict:
    """Delta logLik and Delta BIC between K=3 and K=2 (plus all raw values)."""
    out = {}
    for K, fit in fits.items():
        out[f"logLik_K{K}"] = fit.loglik
        out[f"BIC_K{K}"] = fit.bic(T)
        out[f"ICL_K{K}"] = fit.icl(T)
    if 2 in fits and 3 in fits:
        out["delta_logLik_3_minus_2"] = fits[3].loglik - fits[2].loglik
        out["delta_BIC_3_minus_2"] = fits[3].bic(T) - fits[2].bic(T)
    return out


def mixture_fit_summary_rows(fit: MixtureHMMFit, T: int) -> list[dict]:
    """One row per state-component of a mixture-emission HMM."""
    occ = fit.gamma.mean(axis=0)
    rows = []
    for k in range(fit.K):
        for m in range(fit.M):
            rows.append({
                "K": fit.K,
                "M": fit.M,
                "state": k + 1,
                "component": m + 1,
                "component_weight": fit.weights[k, m],
                "mean": fit.means[k, m],
                "sd": fit.sds[k, m],
                "state_occupation": occ[k],
                "logLik": fit.loglik,
                "BIC": fit.bic(T),
            })
    return rows


def one_step_residuals_gaussian(y: np.ndarray, fit: HMMFit) -> np.ndarray:
    """Approximate one-step standardized prediction residuals."""
    pred = fit.pi.copy()
    residuals = np.zeros(len(y))
    for t, yt in enumerate(y):
        mean = float(np.sum(pred * fit.means))
        second = float(np.sum(pred * (fit.sds ** 2 + fit.means ** 2)))
        var = max(second - mean ** 2, 1e-10)
        residuals[t] = (yt - mean) / np.sqrt(var)
        dens = np.exp(-0.5 * ((yt - fit.means) / fit.sds) ** 2)
        dens /= np.sqrt(2.0 * np.pi) * fit.sds
        filt = pred * np.maximum(dens, 1e-300)
        filt /= filt.sum()
        pred = filt @ fit.A
    return residuals


def residual_acf(residuals: np.ndarray, max_lag: int = 20) -> pd.DataFrame:
    """Sample autocorrelation function for lags 1..max_lag."""
    x = residuals - residuals.mean()
    denom = float(np.dot(x, x))
    rows = []
    for lag in range(1, max_lag + 1):
        acf = float(np.dot(x[:-lag], x[lag:]) / denom) if denom > 0 else np.nan
        rows.append({"lag": lag, "acf": acf})
    return pd.DataFrame(rows)
