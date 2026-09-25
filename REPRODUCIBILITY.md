# Reproducing the analyses

This guide describes the data, fitting procedures and numerical checks for
*Extra Hidden States under Emission Misspecification in Hidden Markov Models*.

The main application analyzes chromosomal allelic imbalance in HCC1143.
Start with [genomics/README.md](genomics/README.md) for the processed inputs,
fixed-parameter evaluation, complete refit and raw microarray preprocessing.
The supplementary elk analyses examine emission sensitivity and temporal
prediction. Their reproduction sequence follows below.

The legacy renv lock concerns the earlier simulation and movement analyses.
For genomic commands, disable automatic activation in the terminal with
`export RENV_CONFIG_AUTOLOADER_ENABLED=false` and follow the genomic dependency
instructions. For the retained legacy environment, restore that variable to
its default before running `renv::restore()` in a separate copy.

## Compile the article and verify package integrity

The article source, bibliography style and figures are self-contained in
`article/`; the supporting-information source is in `supplement/`. With
TeX Live installed, run `latexmk -pdf main.tex` from `article/` and
`latexmk -pdf supporting_information.tex` from `supplement/`.
The SHA-256 manifest lists every distributed file except the manifest itself.
The commands below reproduce the simulation study, genomic analyses and
supplementary movement analyses, with checks of the saved estimates.

## Added four-state genomic comparison

`Rscript genomics/scripts/genomic_g4.R` reproduces the added G4 fit from
40 random starts and nine duplicated-state G3 starts. Checkpoints record
input and reference hashes and are resumed only when those hashes match.
The complete genomic fitting script calls this extension automatically.
`Rscript genomics/scripts/test_genomic_g4.R` cross-checks the G4 likelihood
and the seven unchanged model hashes. Fixed-parameter evaluation scores
are recomputed using R emission densities and the same compiled forward
recursion as the evaluation pipeline; this is not an independent filtering
implementation.
After evaluation and diagnostics, regenerate all genomic tables with
`Rscript numerics/scripts/render_publication_tables.R`.

G4 was added after the original evaluation. Its scores on excluded
chromosomes are retrospective comparisons, as stated in Supporting
Information S3.1. The original model freeze is not extended retroactively.

## Gaussian movement benchmark

Use R 4.5.0 or a compatible version, a compiler supported by Rcpp, and the
package versions in `renv.lock`. From this extracted directory:

```sh
Rscript -e 'renv::restore()'
Rscript numerics/scripts/run_revision_validation.R
Rscript numerics/scripts/refine_revision_sensitivity.R
Rscript numerics/scripts/run_revision_validation.R
Rscript numerics/scripts/check_elk_context.R
Rscript numerics/scripts/render_revision_outputs.R
Rscript numerics/scripts/test_revision.R
```

The first command running the analysis reads existing cached fits. To rerun
from seeds, use a separate extracted copy and remove its generated
`empirical/revision2/fits` directory before executing the sequence. The second
validation call refreshes summaries after the sensitivity refinements.
No training-prefix fit is initialized with a full-series estimate.

The Gaussian benchmark comprises 100 fits in the primary/sensitivity/temporal
validation grid plus five fits to the published movement rates for elk 163.
Each uses 200 initializations, screening, refinement and a convergence audit.
The forward-backward acceleration is checked against the original R recursion
and direct path enumeration. Seeds, settings and session information are saved.

## Results and their scope

`empirical/revision2/` contains the Gaussian benchmark:

- `model_summary.csv`: likelihood, BIC, modes, boundaries and convergence;
- `start_audit.csv`: screening, refined and warm-start solutions;
- `heldout_predictions.csv`: every reserved-observation score and PIT;
- `heldout_summary.csv` and `heldout_by_fold.csv`: predictive summaries;
- `parameters.csv`: component and state parameters;
- `context_summary.csv` and `context_provenance.json`: supplementary-variable checks;
- `decoding_sensitivity.csv`: path stability under alternative variance bounds;
- `published_rate_sensitivity.csv`: the original-rate check;
- `fits/`: all fitted objects, including likelihood traces and start records.

The Gaussian result for elk 163 is qualified by the additional positive-family
comparisons: an ordinary two-state model can have lower BIC. The worked joint
example is elk 287, with heterogeneous short and long movements in a mobile
phase followed by localized movement. This is a retrospective description,
not independently validated behaviour or a general forecasting advantage.

## Joint movement models and emission alternatives

The models are described in Supporting Information S5; the complete
reproduction sequence is in S6. Run these commands
from the package root after the Gaussian benchmark:

```sh
Rscript numerics/scripts/test_movement_review.R
Rscript numerics/scripts/run_movement_review.R
Rscript numerics/scripts/run_movement_review.R joint
Rscript numerics/scripts/run_movement_validation.R
Rscript numerics/scripts/run_movement_validation.R joint
Rscript numerics/scripts/refine_movement_review.R
Rscript numerics/scripts/run_movement_sensitivity.R
Rscript numerics/scripts/refine_movement_starts.R
Rscript numerics/scripts/refine_movement_validation.R
Rscript numerics/scripts/refine_movement_review.R
Rscript numerics/scripts/run_movement_review.R
Rscript numerics/scripts/run_movement_review.R joint
Rscript numerics/scripts/run_movement_validation.R
Rscript numerics/scripts/run_movement_validation.R joint
Rscript numerics/scripts/run_movement_sensitivity.R
Rscript numerics/scripts/run_water_comparison.R
Rscript numerics/scripts/run_extended_gaussian_validation.R
Rscript numerics/scripts/describe_movement_review.R
Rscript numerics/scripts/check_movement_phase.R
Rscript numerics/scripts/test_movement_results.R
Rscript numerics/scripts/render_movement_review.R
```

`empirical/application_review/` contains 160 full-series and 480 training-prefix
fits, 120 resolution/angular sensitivity fits, 48 water-covariate fits, 60
additional Gaussian training fits, and three transition-boundary fits. The
count refers to saved model objects, not the number of optimizer starts.
Use a separate copy and remove this generated directory for a complete seeded
rerun. The downloaded input is checked against its recorded SHA-256 digest.

Primary outputs include `joint_summary.csv`, `rate_summary.csv`,
`joint_validation_blocks.csv`, `joint_validation_with_pit.csv`,
`sensitivity_summary.csv`, `water_comparison.csv`, `elk287_phase_context.csv`,
`component_parameters.csv`, `optimization_audit.csv`, and
`predictive_dependence.csv`. Start-level histories are retained in each fit.
`test_movement_review.R` checks the likelihood independently on small examples;
`test_movement_results.R` rechecks all primary fits and forecasts against the
saved observations, verifies likelihood nesting, and checks the water models.

Optimization codes and projected gradients are separate diagnostics. One
training fit retains a solver warning despite a small projected gradient;
this is documented in the audit, not relabelled as code-zero convergence.
BIC, predictive scores and biological interpretation answer different questions.
The scores are retrospective comparisons; no untouched confirmation sample
is claimed. Sensitivity to rounding resolution is not a GPS-error model.

## Simulation study and numerical optimization checks

The study uses 200 replications at each of 500, 1500 and 5000 observations,
with master seed 20260610. All 600 replications have been rerun using the
original random-start distributions, data-derived nested starts and a
relative likelihood tolerance of 1e-9. The original parameter-generating
values are not used as starts. Reproduce the corrected study with:

```sh
Rscript numerics/scripts/test_fast_em.R
Rscript numerics/scripts/run_optimization_audit.R
Rscript numerics/scripts/refine_reference_simulation.R
Rscript numerics/scripts/test_optimization_audit.R
```

The simulation script resumes saved per-replication checkpoints. To rerun
from scratch, move `numerics/output/optimization_audit/fit_*.rds` to a separate
backup directory first. `original_replications.csv` preserves the original
results; it is not used as the source of the revised tables. The historical
`run_review_simulation.R` implements the original, less thorough fitting
protocol and will overwrite corrected outputs if run directly.

The equal-component genomic comparison additionally has a model-conditional
parametric bootstrap. Its seeds, fits, convergence checks and tail calculation
are saved under `genomics/results/refinement_bootstrap`. Reproduce with:

```sh
Rscript genomics/scripts/genomic_refinement_bootstrap.R
Rscript numerics/scripts/render_publication_tables.R
```

It resumes the 199 saved replicates. Move `observed.rds` and `replicate_*.rds`
to a backup directory to refit them. The model uses two independent sequences
with a common freely estimated initial distribution. This is not a test of
biological state count, and nonrejection does not establish equivalence.

## Earlier movement computations

`numerics/scripts/run_elk_application.R`, `empirical/output/` and the older
empirical figures preserve the earlier Gaussian analysis. The current
supplementary tables are generated by the benchmark and joint-model
workflows above.
Running that older script can overwrite the shared article figure/table
locations; rerun `render_revision_outputs.R` and `render_movement_review.R` afterward to restore the revised
figures and tables. The original EM implementation and original tests are retained.
