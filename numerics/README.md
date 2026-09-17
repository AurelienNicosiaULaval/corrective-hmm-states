# Simulation and supplementary numerical workflows

This directory contains the simulation study and supplementary movement
analyses for *Extra Hidden States under Emission Misspecification in Hidden
Markov Models*. For the main genomic application, see
[genomics/README.md](../genomics/README.md).

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

From the package root:

```sh
Rscript numerics/scripts/test_fast_em.R
Rscript numerics/scripts/run_optimization_audit.R
Rscript numerics/scripts/refine_reference_simulation.R
Rscript numerics/scripts/test_optimization_audit.R
Rscript numerics/tests/testthat.R
```

The corrected simulation uses the same data seeds as the original study,
10 random starts per model, the three best screened candidates, and additional
starts derived from fitted G2, M21 and M22 models. Final fits have relative
log-likelihood gain at most 1e-9, with exact-refinement nesting checks.
The compiled EM is tested against the original R updates. Per-fit parameters,
traces and the original results are saved in `output/optimization_audit`.
The historical `run_review_simulation.R` is retained to document the original
protocol, but running it directly overwrites corrected outputs.

## Main outputs

The simulation scripts write:

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
