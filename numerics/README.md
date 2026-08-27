# Numerical study: corrective hidden states in misspecified HMMs

This directory contains the reproducible numerical study for manuscript 4416781,
*Extra Hidden States under Emission Misspecification in Hidden Markov Models*.

## Idea

The data-generating model is a persistent two-regime aggregate HMM. The first
aggregate regime has a two-component Gaussian mixture emission, while the
second aggregate regime has one Gaussian emission. When fitted with simple
Gaussian-emission HMMs, the heterogeneous first regime can be split into two
fitted states. When fitted with a two-state HMM whose emissions are Gaussian
mixtures, the same correction can be absorbed inside the emission distribution.

## Layout

```text
numerics/
  R/
    simulate.R              true model and simulation
    hmm_em.R                Gaussian-HMM EM
    mixture_hmm_em.R        mixture-emission HMM EM
    diagnostics.R           summaries, BIC, crosstabs and residuals
    plotting.R              ggplot2 article and supplement figures
  scripts/
    run_review_simulation.R
    run_elk_application.R
  tests/
    testthat.R
    testthat/
  output/
```

## Requirements

The HMM implementations use base R matrix operations and do not rely on an
external HMM fitting library. Tables and figures use standard R packages,
including `dplyr`, `ggplot2`, and `patchwork`. The versioned `renv.lock` file
records the complete package environment. Restore it from the repository root
with:

```sh
Rscript -e 'renv::restore()'
```

The elk analysis reads `moveHMM::elk_data` directly in memory and checks its
structure and SHA-256 digest before fitting the models. It neither creates nor
redistributes a CSV containing the original coordinates or derived coordinate
endpoints.

## Main reproduction command

From the project root:

```sh
Rscript numerics/scripts/run_review_simulation.R
Rscript numerics/scripts/run_elk_application.R
Rscript numerics/tests/testthat.R
```

All seeds are fixed from `SEED = 20260610`. The default review run uses
`N_REPS = 200` for each `T` in `(500, 1500, 5000)`, 10 random EM starts for
each replicated Gaussian fit, 10 random EM starts for each replicated
mixture-emission fit and deterministic parallel execution. The replicated fits
screen all 10 starts, then refine the three best starts with relative
log-likelihood tolerance `1e-6`.

For quick local checks only, the number of replications can be reduced:

```sh
REVIEW_N_REPS=2 Rscript numerics/scripts/run_review_simulation.R
```

Do not use reduced-replication outputs for submission.

## Main outputs

`run_review_simulation.R` writes:

- `numerics/output/review_gaussian_fit_summary.csv`;
- `numerics/output/review_mixture_fit_summary.csv`;
- `numerics/output/review_model_comparison.csv`;
- `numerics/output/review_replications.csv`;
- `numerics/output/review_summary_by_T.csv`;
- `numerics/output/review_metadata.json`;
- `numerics/output/residual_acf_fit2.csv`;
- `numerics/output/residual_acf_fit3.csv`;
- `article/tables/simulation_table.tex`;
- `article/tables/replication_summary_table.tex`;
- `supplement/tables/review_summary_table.tex`;
- `article/figures/simulation_density.{pdf,png}`;
- `article/figures/state_decoding_heatmap.{pdf,png}`;
- `article/figures/posterior_corrective_state.{pdf,png}`;
- `article/figures/bic_boxplot.{pdf,png}`;
- `article/figures/bic_by_T.{pdf,png}`;
- `article/figures/residual_acf.{pdf,png}`.

## Model comparison

The main comparison uses:

- Gaussian-emission HMM, `K = 2`;
- Gaussian-emission HMM, `K = 3`;
- mixture-emission HMM, `K = 2`, `M = 2` Gaussian components per state.

The BIC parameter count for the mixture-emission HMM is
`(K - 1) + K(K - 1) + K((M - 1) + 2M)`, covering the initial distribution,
transitions, mixture weights, means and standard deviations.

The numerical interpretation should always be read from the generated outputs.
In the reference run, the intended diagnostic pattern is that the simple
Gaussian family favours an additional fitted state, while the enriched
two-state mixture-emission model absorbs the split into the emission model.
