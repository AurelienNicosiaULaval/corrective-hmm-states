# Corrective hidden states in misspecified HMMs

This repository contains the manuscript, supporting information, reproducible
numerical code, and generated outputs for the *Stat* manuscript:

*Extra Hidden States as Kullback-Leibler Corrections under Emission
Misspecification in Hidden Markov Models*.

## Repository contents

```text
article/
  main.tex                  compact LaTeX manuscript
  main.pdf                  compiled manuscript PDF
  references.bib            BibTeX bibliography
  figures/                  generated numerical figures
  tables/                   generated LaTeX tables
supplement/
  supporting_information.tex
  supporting_information.pdf
  tables/                   supporting-information tables
numerics/
  scripts/run_review_simulation.py
  src/                      simulation, HMM EM, diagnostics and figures
  tests/                    numerical tests
  output/                   generated CSV/JSON/LaTeX outputs
tools/
  validate_submission.py
  build_submission_package.py
```

## Manuscript

The manuscript is intentionally compact and centers on Theorem 5.1. Under
sufficient separation conditions, a refined HMM reproduces the observed law
exactly, whereas a non-refined HMM with the aggregate number of states remains
separated from the truth in Kullback-Leibler rate.

The supporting information gives the mixture-emission EM details, numerical
replication information, supplementary diagnostics, and reproducibility notes.

## Numerical reproduction

From the project root:

```bash
python3 numerics/scripts/run_review_simulation.py
python3 -m pytest numerics/tests -q
```

The simulation script writes:

- PDF and PNG figures in `article/figures/`;
- the main manuscript table in `article/tables/simulation_table.tex`;
- the replication table in `article/tables/replication_summary_table.tex`;
- the supporting-information table in `supplement/tables/review_summary_table.tex`;
- CSV and JSON outputs in `numerics/output/`.

No external data are used.

## Compilation

From `article/`:

```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

From `supplement/`:

```bash
pdflatex supporting_information.tex
pdflatex supporting_information.tex
```

## Local validation

From the project root:

```bash
python3 tools/validate_submission.py
python3 tools/build_submission_package.py
```

The validator checks required files, public-facing text markers, 200
replications per sample size, at least 10 EM starts for each replicated fit,
deterministic start screening, and a short numerical smoke test including the
mixture-emission HMM.

## Submission Archive

`tools/build_submission_package.py` rebuilds
`submission/stat_hmm_submission_package/` and
`submission/stat_hmm_submission_package.zip` with the manuscript, supporting
information, LaTeX sources, figures, tables, numerical code, and reproducible
outputs.
