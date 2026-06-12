"""Pre-submission review numerics.

This script extends the original illustration by adding:

- BIC selection frequencies over T in {500, 1500, 5000};
- an enriched two-state HMM with two Gaussian components per state;
- one-step residual ACF diagnostics for the main series;
- article and supporting-information tables.

Run from the project root:
    python3 numerics/scripts/run_review_simulation.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# The replicated EM fits are many small linear algebra problems. Limiting BLAS
# threads avoids oversubscription when several Python workers run in parallel.
for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "numerics"))
os.environ.setdefault("MPLBACKEND", "Agg")

from src.diagnostics import (fit_summary_rows, mixture_fit_summary_rows,
                             one_step_residuals_gaussian, residual_acf)
from src.hmm_em import fit_hmm, validate_fit
from src.mixture_hmm_em import fit_mixture_hmm, validate_mixture_fit
from src.plotting import (plot_bic_boxplot, plot_bic_by_T, plot_decoding_heatmap,
                          plot_density_review, plot_posterior_corrective_state,
                          plot_residual_acf)
from src.simulate import (GAMMA0, MEANS_REFINED, MIX_WEIGHTS_REGIME1,
                          SDS_REFINED, simulate_hmm, stationary_distribution)
from src.diagnostics import crosstab_true_vs_fitted


SEED = 20260610
T_MAIN = 1500
T_GRID = (500, 1500, 5000)
N_REPS = int(os.environ.get("REVIEW_N_REPS", "200"))
N_STARTS_MAIN = int(os.environ.get("REVIEW_N_STARTS_MAIN", "10"))
N_STARTS_REPS = int(os.environ.get("REVIEW_N_STARTS_REPS", "10"))
N_STARTS_MIX_MAIN = int(os.environ.get("REVIEW_N_STARTS_MIX_MAIN", "10"))
N_STARTS_MIX_REPS = int(os.environ.get("REVIEW_N_STARTS_MIX_REPS", "10"))
MAX_ITER_REPS = int(os.environ.get("REVIEW_MAX_ITER_REPS", "100"))
SCREEN_ITER_REPS = int(os.environ.get("REVIEW_SCREEN_ITER_REPS", "25"))
REFINE_TOP_REPS = int(os.environ.get("REVIEW_REFINE_TOP_REPS", "3"))
SCREEN_ITER_MIX_REPS = int(os.environ.get("REVIEW_SCREEN_ITER_MIX_REPS", "30"))
REFINE_TOP_MIX_REPS = int(os.environ.get("REVIEW_REFINE_TOP_MIX_REPS", "3"))
N_WORKERS = int(os.environ.get(
    "REVIEW_N_WORKERS",
    str(max(1, min(8, (os.cpu_count() or 2) - 1)))
))

OUT_DIR = ROOT / "numerics" / "output"
FIG_DIR = ROOT / "article" / "figures"
TAB_DIR = ROOT / "article" / "tables"
SUPP_TAB_DIR = ROOT / "supplement" / "tables"
PARTIAL_REPS_PATH = OUT_DIR / "review_replications.partial.csv"
PARTIAL_META_PATH = OUT_DIR / "review_replications.partial.json"


def main() -> None:
    t_start = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    SUPP_TAB_DIR.mkdir(parents=True, exist_ok=True)

    print("Fitting reference series models", flush=True)
    y, _, refined = simulate_hmm(T_MAIN, SEED)
    fit2 = fit_hmm(y, 2, seed=SEED + 2, n_starts=N_STARTS_MAIN)
    fit3 = fit_hmm(y, 3, seed=SEED + 3, n_starts=N_STARTS_MAIN)
    mix2 = fit_mixture_hmm(y, 2, M=2, seed=SEED + 20,
                           n_starts=N_STARTS_MIX_MAIN)
    validate_fit(fit2, T_MAIN)
    validate_fit(fit3, T_MAIN)
    validate_mixture_fit(mix2, T_MAIN)

    gaussian_summary = pd.DataFrame(
        fit_summary_rows(fit2, T_MAIN) + fit_summary_rows(fit3, T_MAIN))
    gaussian_summary.to_csv(OUT_DIR / "review_gaussian_fit_summary.csv", index=False)
    mixture_summary = pd.DataFrame(mixture_fit_summary_rows(mix2, T_MAIN))
    mixture_summary.to_csv(OUT_DIR / "review_mixture_fit_summary.csv", index=False)

    tab2 = crosstab_true_vs_fitted(refined, fit2, use="viterbi")
    tab3 = crosstab_true_vs_fitted(refined, fit3, use="viterbi")
    tab2.to_csv(OUT_DIR / "crosstab_fit2.csv")
    tab3.to_csv(OUT_DIR / "crosstab_fit3.csv")

    plot_density_review(y, fit2, fit3, mix2, FIG_DIR)
    plot_decoding_heatmap(tab2, tab3, FIG_DIR)
    corrective_idx = int(np.argmin(np.abs(fit3.means - MEANS_REFINED[1])))
    plot_posterior_corrective_state(y, fit3, corrective_idx, FIG_DIR)
    acf2 = residual_acf(one_step_residuals_gaussian(y, fit2))
    acf3 = residual_acf(one_step_residuals_gaussian(y, fit3))
    acf2.to_csv(OUT_DIR / "residual_acf_fit2.csv", index=False)
    acf3.to_csv(OUT_DIR / "residual_acf_fit3.csv", index=False)
    plot_residual_acf(acf2, acf3, FIG_DIR)
    print("Reference series models and figures done", flush=True)

    rep_rows = _run_replications(t_start)

    reps = pd.DataFrame(rep_rows)
    reps = reps.sort_values(["T", "rep"]).reset_index(drop=True)
    reps.to_csv(OUT_DIR / "review_replications.csv", index=False)
    summary = reps.groupby("T", as_index=False).agg(
        n_replications=("rep", "count"),
        mean_delta_BIC_gauss_K3_minus_K2=("delta_BIC_gauss_K3_minus_K2", "mean"),
        prop_gaussian_BIC_favours_K3=("gaussian_BIC_favours_K3", "mean"),
        mean_delta_BIC_mix_K2_minus_gauss_K3=("delta_BIC_mix_K2_minus_gauss_K3", "mean"),
        prop_enriched_BIC_favours_K2_mix=("enriched_BIC_favours_K2_mix", "mean"),
    )
    summary.to_csv(OUT_DIR / "review_summary_by_T.csv", index=False)
    plot_bic_boxplot(
        reps.loc[reps["T"] == T_MAIN, "delta_BIC_gauss_K3_minus_K2"].to_numpy(),
        FIG_DIR)
    plot_bic_by_T(summary, FIG_DIR)

    model_rows = _model_rows(fit2, fit3, mix2, T_MAIN)
    model_table = pd.DataFrame(model_rows)
    model_table.to_csv(OUT_DIR / "review_model_comparison.csv", index=False)
    _write_article_table(model_rows)
    _write_replication_tables(summary)

    metadata = {
        "seed": SEED,
        "T_main": T_MAIN,
        "T_grid": list(T_GRID),
        "n_replications_per_T": N_REPS,
        "n_starts_main": N_STARTS_MAIN,
        "n_starts_replications": N_STARTS_REPS,
        "n_starts_mixture_main": N_STARTS_MIX_MAIN,
        "n_starts_mixture_replications": N_STARTS_MIX_REPS,
        "n_workers": N_WORKERS,
        "max_iter_replications": MAX_ITER_REPS,
        "screen_iter_replications": SCREEN_ITER_REPS,
        "refine_top_replications": REFINE_TOP_REPS,
        "screen_iter_mixture_replications": SCREEN_ITER_MIX_REPS,
        "refine_top_mixture_replications": REFINE_TOP_MIX_REPS,
        "true_transition": GAMMA0.tolist(),
        "true_stationary": stationary_distribution(GAMMA0).tolist(),
        "true_mixture_weights_regime1": MIX_WEIGHTS_REGIME1.tolist(),
        "true_refined_means": MEANS_REFINED.tolist(),
        "true_refined_sds": SDS_REFINED.tolist(),
        "main_model_comparison": model_rows,
        "selection_summary": summary.to_dict(orient="records"),
        "runtime_seconds": round(time.time() - t_start, 1),
    }
    (OUT_DIR / "review_metadata.json").write_text(json.dumps(metadata, indent=2))
    PARTIAL_REPS_PATH.unlink(missing_ok=True)
    PARTIAL_META_PATH.unlink(missing_ok=True)
    print(json.dumps(metadata, indent=2))


def _fit_replication(args: tuple[int, int]) -> dict:
    T, rep = args
    seed_r = SEED + 1_000_000 + 10_000 * T + rep
    y_r, _, _ = simulate_hmm(T, seed_r)
    f2 = fit_hmm(y_r, 2, seed=seed_r + 2, n_starts=N_STARTS_REPS,
                 max_iter=MAX_ITER_REPS, screen_iter=SCREEN_ITER_REPS,
                 refine_top=REFINE_TOP_REPS)
    f3 = fit_hmm(y_r, 3, seed=seed_r + 3, n_starts=N_STARTS_REPS,
                 max_iter=MAX_ITER_REPS, screen_iter=SCREEN_ITER_REPS,
                 refine_top=REFINE_TOP_REPS)
    fm = fit_mixture_hmm(y_r, 2, M=2, seed=seed_r + 20,
                         n_starts=N_STARTS_MIX_REPS,
                         max_iter=MAX_ITER_REPS,
                         screen_iter=SCREEN_ITER_MIX_REPS,
                         refine_top=REFINE_TOP_MIX_REPS)
    validate_fit(f2, T)
    validate_fit(f3, T)
    validate_mixture_fit(fm, T)
    return {
        "T": T,
        "rep": rep,
        "seed": seed_r,
        "logLik_gauss_K2": f2.loglik,
        "logLik_gauss_K3": f3.loglik,
        "logLik_mix_K2_M2": fm.loglik,
        "BIC_gauss_K2": f2.bic(T),
        "BIC_gauss_K3": f3.bic(T),
        "BIC_mix_K2_M2": fm.bic(T),
        "delta_BIC_gauss_K3_minus_K2": f3.bic(T) - f2.bic(T),
        "delta_BIC_mix_K2_minus_gauss_K3": fm.bic(T) - f3.bic(T),
        "gaussian_BIC_favours_K3": f3.bic(T) < f2.bic(T),
        "enriched_BIC_favours_K2_mix": fm.bic(T) < f3.bic(T),
        "mix_state1_mean": fm.state_means()[0],
        "mix_state2_mean": fm.state_means()[1],
        "mix_state1_occ": fm.gamma.mean(axis=0)[0],
        "mix_state2_occ": fm.gamma.mean(axis=0)[1],
    }


def _run_replications(t_start: float) -> list[dict]:
    tasks_all = [(T, rep) for T in sorted(T_GRID, reverse=True) for rep in range(N_REPS)]
    rows = _load_partial_rows()
    completed = {(int(row["T"]), int(row["rep"])) for row in rows}
    tasks = [task for task in tasks_all if task not in completed]
    total = len(tasks_all)
    done_by_T = {T: 0 for T in T_GRID}
    for T, _ in completed:
        if T in done_by_T:
            done_by_T[T] += 1
    if rows:
        print(f"Resuming from {len(rows)}/{total} completed replications",
              flush=True)
    _write_partial_rows(rows)
    print(f"Running {len(tasks)} remaining replications out of {total} "
          f"with {N_WORKERS} worker(s)", flush=True)
    if not tasks:
        return rows

    if N_WORKERS == 1:
        for task in tasks:
            row = _fit_replication(task)
            rows.append(row)
            done_by_T[row["T"]] += 1
            _write_partial_rows(rows)
            if len(rows) % 25 == 0 or len(rows) == total:
                progress = ", ".join(
                    f"T={T}: {done_by_T[T]}/{N_REPS}" for T in T_GRID)
                print(f"{len(rows)}/{total} replications done [{progress}] "
                      f"({time.time() - t_start:.0f}s elapsed)", flush=True)
        return rows

    with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        future_to_task = {executor.submit(_fit_replication, task): task for task in tasks}
        for future in as_completed(future_to_task):
            row = future.result()
            rows.append(row)
            done_by_T[row["T"]] += 1
            _write_partial_rows(rows)
            if len(rows) % 25 == 0 or len(rows) == total:
                progress = ", ".join(
                    f"T={T}: {done_by_T[T]}/{N_REPS}" for T in T_GRID)
                print(f"{len(rows)}/{total} replications done [{progress}] "
                      f"({time.time() - t_start:.0f}s elapsed)", flush=True)
    return rows


def _partial_config() -> dict:
    return {
        "seed": SEED,
        "T_grid": list(T_GRID),
        "n_replications_per_T": N_REPS,
        "n_starts_replications": N_STARTS_REPS,
        "n_starts_mixture_replications": N_STARTS_MIX_REPS,
        "max_iter_replications": MAX_ITER_REPS,
        "screen_iter_replications": SCREEN_ITER_REPS,
        "refine_top_replications": REFINE_TOP_REPS,
        "screen_iter_mixture_replications": SCREEN_ITER_MIX_REPS,
        "refine_top_mixture_replications": REFINE_TOP_MIX_REPS,
    }


def _load_partial_rows() -> list[dict]:
    config = _partial_config()
    if not PARTIAL_REPS_PATH.exists() or not PARTIAL_META_PATH.exists():
        PARTIAL_META_PATH.write_text(json.dumps(config, indent=2))
        return []
    existing_config = json.loads(PARTIAL_META_PATH.read_text())
    if existing_config != config:
        print("Ignoring partial replication file because its configuration differs",
              flush=True)
        PARTIAL_REPS_PATH.unlink(missing_ok=True)
        PARTIAL_META_PATH.write_text(json.dumps(config, indent=2))
        return []
    partial = pd.read_csv(PARTIAL_REPS_PATH)
    partial = partial.loc[
        partial["T"].isin(T_GRID)
        & partial["rep"].between(0, N_REPS - 1)
    ].drop_duplicates(["T", "rep"], keep="last")
    partial = partial.sort_values(["T", "rep"]).reset_index(drop=True)
    return partial.to_dict(orient="records")


def _write_partial_rows(rows: list[dict]) -> None:
    PARTIAL_META_PATH.write_text(json.dumps(_partial_config(), indent=2))
    if not rows:
        return
    partial = pd.DataFrame(rows)
    partial = partial.sort_values(["T", "rep"]).drop_duplicates(
        ["T", "rep"], keep="last")
    tmp_path = PARTIAL_REPS_PATH.with_suffix(".partial.csv.tmp")
    partial.to_csv(tmp_path, index=False)
    tmp_path.replace(PARTIAL_REPS_PATH)


def _model_rows(fit2, fit3, mix2, T):
    rows = [
        {
            "model": "Gaussian HMM",
            "emission": "one Gaussian per state",
            "K": 2,
            "logLik": fit2.loglik,
            "BIC": fit2.bic(T),
            "diagnostic": "merges the heterogeneous regime",
        },
        {
            "model": "Gaussian HMM",
            "emission": "one Gaussian per state",
            "K": 3,
            "logLik": fit3.loglik,
            "BIC": fit3.bic(T),
            "diagnostic": "splits the heterogeneous regime",
        },
        {
            "model": "Mixture-emission HMM",
            "emission": "two Gaussians per state",
            "K": 2,
            "logLik": mix2.loglik,
            "BIC": mix2.bic(T),
            "diagnostic": "absorbs the split in the emission",
        },
    ]
    return rows


def _write_article_table(rows) -> None:
    specification_label = {
        ("Gaussian HMM", "one Gaussian per state"): "Gaussian, one Gaussian/state",
        ("Mixture-emission HMM", "two Gaussians per state"): "Mixture, two Gaussians/state",
    }
    diagnostic_label = {
        "merges the heterogeneous regime": "merged regime",
        "splits the heterogeneous regime": "split regime",
        "absorbs the split in the emission": "emission split",
    }
    lines = [
        r"% Generated by numerics/scripts/run_review_simulation.py.",
        r"\begin{table}[t]",
        r"\centering",
        r"\small",
        r"\caption{Main numerical diagnostic, $T=1500$. Lower BIC is better.}",
        r"\label{tab:sim}",
        r"\begin{tabular}{@{}lrrrp{0.24\linewidth}@{}}",
        r"\toprule",
        r"Specification & $K$ & logLik & BIC & interpretation\\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            f"{specification_label.get((row['model'], row['emission']), row['model'])} & "
            f"{row['K']} & {row['logLik']:.1f} & {row['BIC']:.1f} & "
            f"{diagnostic_label.get(row['diagnostic'], row['diagnostic'])}\\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    (TAB_DIR / "simulation_table.tex").write_text("\n".join(lines))


def _write_replication_tables(summary: pd.DataFrame) -> None:
    lines = [
        r"% Generated by numerics/scripts/run_review_simulation.py.",
        r"\begin{table}[t]",
        r"\centering",
        r"\small",
        r"\caption{BIC selection frequencies over fixed-seed replications. Lower BIC is better; $\Delta_{\rm G}=\mathrm{BIC}_{G,3}-\mathrm{BIC}_{G,2}$ and $\Delta_{\rm M}=\mathrm{BIC}_{M,2}-\mathrm{BIC}_{G,3}$.}",
        r"\label{tab:replications}",
        r"\begin{tabular}{rrrrr}",
        r"\toprule",
        r"$T$ & reps & mean $\Delta_{\rm G}$ & Pr$(\Delta_{\rm G}<0)$ & Pr$(\Delta_{\rm M}<0)$\\",
        r"\midrule",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"{int(row['T'])} & {int(row['n_replications'])} & "
            f"{row['mean_delta_BIC_gauss_K3_minus_K2']:.1f} & "
            f"{row['prop_gaussian_BIC_favours_K3']:.3f} & "
            f"{row['prop_enriched_BIC_favours_K2_mix']:.3f}\\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    table_text = "\n".join(lines)
    (TAB_DIR / "replication_summary_table.tex").write_text(table_text)
    (SUPP_TAB_DIR / "review_summary_table.tex").write_text(table_text)


if __name__ == "__main__":
    main()
