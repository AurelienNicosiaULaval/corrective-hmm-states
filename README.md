# Extra Hidden States under Emission Misspecification in Hidden Markov Models

This repository contains the revised manuscript, supporting information,
analysis code, fixed seeds, generated outputs, and validation tests for
manuscript 4416781, currently under revision for *Statistica Neerlandica*.

The numerical work has two parts: a fixed-seed simulation study and an
empirical illustration using the public `elk_data` object distributed with the
R package `moveHMM`. The complete computational workflow is implemented in R,
and all figures are generated with `ggplot2`.

## Repository contents

```text
article/                         revised manuscript, figures and tables
supplement/                      proofs, computational details and diagnostics
numerics/R/                      HMM implementations and diagnostics
numerics/scripts/                simulation and elk analysis drivers
numerics/tests/                  testthat validation suite
numerics/output/                 fixed-seed simulation outputs
empirical/output/                derived elk analysis outputs
empirical/figures/               elk diagnostic figures
empirical/DATA_SOURCE.md         data provenance and object digest
renv.lock                        locked R package environment
```

## Requirements

The archived analyses were produced with R 4.5.0 and `moveHMM` 1.10. Package
versions are recorded in `renv.lock`. From the repository root, restore the
environment with:

```sh
Rscript -e 'renv::restore()'
```

A LaTeX distribution is required only to rebuild the document PDFs.

## Reproduce the analyses

From the repository root:

```sh
Rscript numerics/scripts/run_review_simulation.R
Rscript numerics/scripts/run_elk_application.R
Rscript numerics/tests/testthat.R
Rscript tools/validate_repository.R
```

The full simulation uses 200 replications at each sample size. A reduced
implementation check can be run with:

```sh
REVIEW_N_REPS=2 Rscript numerics/scripts/run_review_simulation.R
```

Reduced outputs must not replace the archived 200-replication results.

## Data provenance

The original elk coordinate table is not redistributed. The empirical script
loads `moveHMM::elk_data` directly into memory and verifies its structure and
SHA-256 object digest before fitting any model. No CSV containing the original
coordinates or derived coordinate endpoints is created. The CSV files in
`empirical/output/` contain derived analysis results only. Full provenance and
source citations are recorded in `empirical/DATA_SOURCE.md`.

## Compile the documents

From `article/`:

```sh
latexmk -pdf main.tex
```

From `supplement/`:

```sh
latexmk -pdf supporting_information.tex
```

Precompiled PDFs are included for review convenience. Numerical interpretation
should be based on the archived CSV and JSON outputs.
