"""Figures for the article. All figures are saved in PDF and high-res PNG."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .hmm_em import HMMFit
from .mixture_hmm_em import MixtureHMMFit

DPI = 300


def _save(fig, fig_dir: Path, name: str) -> None:
    fig.tight_layout()
    fig.savefig(fig_dir / f"{name}.pdf")
    fig.savefig(fig_dir / f"{name}.png", dpi=DPI)
    plt.close(fig)


def _norm_pdf(x, mu, sd):
    return np.exp(-0.5 * ((x - mu) / sd) ** 2) / (np.sqrt(2.0 * np.pi) * sd)


def plot_density(y: np.ndarray, fits: dict[int, HMMFit], fig_dir: Path) -> None:
    """Histogram of observations with fitted marginal mixture densities."""
    grid = np.linspace(y.min() - 0.8, y.max() + 0.8, 600)
    fig, ax = plt.subplots(figsize=(6.2, 3.7))
    ax.hist(y, bins=60, density=True, alpha=0.35, color="grey",
            label="observations")
    styles = {2: ("--", "C1"), 3: ("-", "C0"), 4: (":", "C2")}
    for K in sorted(fits):
        fit = fits[K]
        occ = fit.gamma.mean(axis=0)
        dens = sum(occ[j] * _norm_pdf(grid, fit.means[j], fit.sds[j])
                   for j in range(K))
        ls, c = styles.get(K, ("-", None))
        ax.plot(grid, dens, ls, color=c, lw=2, label=f"Gaussian HMM, K={K}")
    ax.set_xlabel("Observation")
    ax.set_ylabel("Density")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, fig_dir, "simulation_density")


def plot_density_review(y: np.ndarray, fit2: HMMFit, fit3: HMMFit,
                        mix2: MixtureHMMFit, fig_dir: Path) -> None:
    """Histogram with Gaussian HMM and enriched-emission HMM fitted marginals."""
    grid = np.linspace(y.min() - 0.8, y.max() + 0.8, 600)
    fig, ax = plt.subplots(figsize=(6.2, 3.7))
    ax.hist(y, bins=60, density=True, alpha=0.32, color="grey",
            label="observations")
    for fit, label, ls, color in [
        (fit2, "Gaussian HMM, K=2", "--", "C1"),
        (fit3, "Gaussian HMM, K=3", "-", "C0"),
    ]:
        occ = fit.gamma.mean(axis=0)
        dens = sum(occ[j] * _norm_pdf(grid, fit.means[j], fit.sds[j])
                   for j in range(fit.K))
        ax.plot(grid, dens, ls, color=color, lw=2, label=label)
    occ = mix2.gamma.mean(axis=0)
    dens = np.zeros_like(grid)
    for k in range(mix2.K):
        for m in range(mix2.M):
            dens += occ[k] * mix2.weights[k, m] * _norm_pdf(
                grid, mix2.means[k, m], mix2.sds[k, m])
    ax.plot(grid, dens, ":", color="C2", lw=2.2,
            label="Mixture-emission HMM, K=2")
    ax.set_xlabel("Observation")
    ax.set_ylabel("Density")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, fig_dir, "simulation_density")


def plot_decoding_heatmap(tab2: pd.DataFrame, tab3: pd.DataFrame,
                          fig_dir: Path) -> None:
    """Row-proportion heatmaps: true refined state vs fitted state, K=2 and K=3."""
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for ax, tab, K in [(axes[0], tab2, 2), (axes[1], tab3, 3)]:
        prop = tab.div(tab.sum(axis=1), axis=0)
        im = ax.imshow(prop.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(prop.shape[1]), [str(c) for c in prop.columns])
        ax.set_yticks(range(prop.shape[0]), list(prop.index))
        ax.set_xlabel(f"Fitted state (K={K})")
        if K == 2:
            ax.set_ylabel("True refined state")
        for i in range(prop.shape[0]):
            for j in range(prop.shape[1]):
                v = prop.values[i, j]
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=8, color="white" if v > 0.5 else "black")
    fig.colorbar(im, ax=axes, fraction=0.035, pad=0.03)
    fig.savefig(fig_dir / "state_decoding_heatmap.pdf", bbox_inches="tight")
    fig.savefig(fig_dir / "state_decoding_heatmap.png", dpi=DPI,
                bbox_inches="tight")
    plt.close(fig)


def plot_posterior_corrective_state(y: np.ndarray, fit3: HMMFit,
                                    corrective_idx: int, fig_dir: Path) -> None:
    """Posterior probability of the corrective state as a function of Y_t."""
    fig, ax = plt.subplots(figsize=(6.2, 3.5))
    ax.scatter(y, fit3.gamma[:, corrective_idx], s=6, alpha=0.30,
               color="C0", edgecolors="none")
    # Binned mean for readability
    bins = np.quantile(y, np.linspace(0, 1, 41))
    idx = np.clip(np.digitize(y, bins) - 1, 0, 39)
    centers, means_ = [], []
    for b in range(40):
        m = idx == b
        if m.sum() >= 5:
            centers.append(y[m].mean())
            means_.append(fit3.gamma[m, corrective_idx].mean())
    ax.plot(centers, means_, color="C3", lw=2, label="binned mean")
    ax.set_xlabel(r"Observation $Y_t$")
    ax.set_ylabel("Posterior prob. of corrective state")
    ax.set_ylim(-0.03, 1.03)
    ax.legend(frameon=False, fontsize=8)
    _save(fig, fig_dir, "posterior_corrective_state")


def plot_bic_boxplot(delta_bic: np.ndarray, fig_dir: Path) -> None:
    """Boxplot of BIC(K=3) - BIC(K=2) across replications."""
    fig, ax = plt.subplots(figsize=(4.4, 3.5))
    ax.boxplot(delta_bic, vert=True, widths=0.4)
    ax.axhline(0.0, color="C3", lw=1, ls="--")
    ax.set_xticks([1], ["replications"])
    ax.set_ylabel(r"BIC$(K{=}3)$ $-$ BIC$(K{=}2)$")
    prop = float(np.mean(delta_bic < 0))
    ax.set_title(f"BIC favours K=3 in {100 * prop:.0f}% of replications",
                 fontsize=9)
    _save(fig, fig_dir, "bic_boxplot")


def plot_bic_by_T(summary: pd.DataFrame, fig_dir: Path) -> None:
    """Selection frequencies by sample size for Gaussian and enriched models."""
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    width = 0.35
    x = np.arange(len(summary))
    ax.bar(x - width / 2, summary["prop_gaussian_BIC_favours_K3"],
           width=width, color="C0", label="Gaussian: BIC favours K=3")
    ax.bar(x + width / 2, summary["prop_enriched_BIC_favours_K2_mix"],
           width=width, color="C2", label="Enriched: BIC favours K=2 mixture")
    ax.set_xticks(x, [str(int(v)) for v in summary["T"]])
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("T")
    ax.set_ylabel("Selection frequency")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, fig_dir, "bic_by_T")


def plot_residual_acf(acf2: pd.DataFrame, acf3: pd.DataFrame, fig_dir: Path) -> None:
    """ACF of one-step standardized prediction residuals."""
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    ax.axhline(0.0, color="black", lw=0.8)
    ax.plot(acf2["lag"], acf2["acf"], marker="o", ms=3, lw=1.2,
            color="C1", label="Gaussian HMM, K=2")
    ax.plot(acf3["lag"], acf3["acf"], marker="o", ms=3, lw=1.2,
            color="C0", label="Gaussian HMM, K=3")
    ax.set_xlabel("Lag")
    ax.set_ylabel("Residual ACF")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, fig_dir, "residual_acf")
