# Extra Hidden States under Emission Misspecification in Hidden Markov Models

Data, code and computational results accompanying
*Extra Hidden States under Emission Misspecification in Hidden Markov Models*.

The main application examines chromosomal allelic imbalance in the HCC1143
breast cancer cell line and its matched normal cell line, HCC1143BL.
Additional analyses of elk movement, including sensitivity analyses and
negative results, are presented in Supporting Information S4–S5.

The repository includes the manuscript and Supporting Information, processed
genomic data, analysis scripts, fitted models, 600 simulation replicates and
199 parametric-bootstrap replicates.

## Start here

| Material | Location |
| --- | --- |
| Current manuscript | [PDF](article/main.pdf), [LaTeX source](article/main.tex) |
| Supporting Information | [PDF](supplement/supporting_information.pdf), [LaTeX source](supplement/supporting_information.tex) |
| Genomic data and analysis guide | [genomics/README.md](genomics/README.md) |
| Full reproduction instructions | [REPRODUCIBILITY.md](REPRODUCIBILITY.md) |
| Data download URLs and SHA-256 hashes | [genomics/source_manifest.json](genomics/source_manifest.json) |
| Version history | [CHANGELOG.md](CHANGELOG.md) |
| Downloadable materials | [Releases](https://github.com/AurelienNicosiaULaval/corrective-hmm-states/releases/latest) |

The scientific question is whether extra Gaussian HMM states describe
additional persistent regimes or compensate for an overly simple emission
distribution. In the genomic example, an unphased allele fraction can alternate
between low and high values within the same regional imbalance. This gives a
measurement-based reason for paired emission components.

## Where are the data?

The original SNP-array measurements are public in
[NCBI GEO, GSE13372](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13372),
from [Chiang et al. (2009), Nature Methods, 6, 99–103](https://doi.org/10.1038/nmeth.1276).
The selected arrays are GSM337641/GSM337662 (tumor/normal) and
GSM337642/GSM337663 (technical repeat tumor/normal).

The processed inputs used by the HMMs are included here. Each row is one
retained SNP marker with its genomic position and allele-fraction measurement.
Coordinates use the hg18 reference genome. These are repeated measurements
of one cell-line pair, not independent patients.

| Processed input | Use | Retained markers |
| --- | --- | ---: |
| [development_chr1.csv](genomics/data/genomics/development_chr1.csv) | Chromosome 1: model fitting | 20,697 |
| [validation_chr12_17.csv](genomics/data/genomics/validation_chr12_17.csv) | Chromosomes 12–17: held-out evaluation | 52,294 |
| [confirmation_chr18_22.csv](genomics/data/genomics/confirmation_chr18_22.csv) | Chromosomes 18–22: later held-out evaluation | 25,413 |
| [replication_chr1.csv](genomics/data/genomics/replication_chr1.csv) | Chromosome 1: technical repeat | 20,937 |

The `confirmation` filename identifies a later held-out block in an exploratory
analysis, not a preregistered confirmatory study. Published sequencing segment
estimates used as an independent assay are also supplied in
[Chiang2009_sequencing_segments_hg18.csv](genomics/data/genomics/Chiang2009_sequencing_segments_hg18.csv).
Their original source is the supplementary material of Chiang et al. (2009).
Raw CEL files are downloaded from GEO when regenerating the preprocessing;
they are not duplicated in this repository. Download locations and hashes,
processing details, source references and dependency versions are documented
in the [genomic guide](genomics/README.md).

For the supplementary elk analyses, see [empirical/DATA_SOURCE.md](empirical/DATA_SOURCE.md),
[the contextual-data provenance](empirical/revision2/context_provenance.json)
and Supporting Information S4–S5. The original coordinates are loaded from
`moveHMM::elk_data`; the retained derived results do not replace that source.

## Repository structure

```text
article/                         current genomic manuscript, figures and tables
supplement/                      proofs, genomics and supplementary elk analyses
genomics/data/genomics/           processed genomic inputs and signal exports
genomics/scripts/                 genomic models, evaluation, preprocessing and tests
genomics/results/                 fitted models, diagnostics and bootstrap replicates
genomics/audit/                   computational provenance and package inventory
numerics/R/                      HMM and mixture-emission EM implementations
numerics/output/optimization_audit/ 600 corrected simulation checkpoints
numerics/scripts/                simulation and supplementary movement workflows
empirical/                       retained supplementary movement results
tools/validate_repository.R      integrity, data-provenance and completeness checks
manifest_sha256.json             SHA-256 digest of every distributed package file
```

## Verify the distributed results

Run commands from the repository root. R 4.5.0 and an Rcpp-compatible C++
compiler were used. For the checks below, install `Rcpp`, `data.table`,
`testthat`, `jsonlite`, `digest` and `moveHMM` (the archived movement data use
version 1.10). Disable the legacy project autoloader before genomic commands:

```sh
export RENV_CONFIG_AUTOLOADER_ENABLED=false
Rscript tools/validate_repository.R
Rscript numerics/tests/testthat.R
Rscript genomics/scripts/test_genomic_results.R
Rscript genomics/scripts/test_symmetric_genomic_hmm.R
Rscript numerics/scripts/test_optimization_audit.R
```

These checks use saved inputs and results; they do not refit the complete
simulation, genomic or movement analyses. The genomic tests independently
check likelihood calculations and the stored model summaries. The optimization
checks examine all 600 simulation checkpoints and 199 bootstrap replicates.
The repository workflow runs these checks on GitHub.

For fixed-parameter evaluation, complete refits, figure generation and raw
array preprocessing, follow [genomics/README.md](genomics/README.md) and
[REPRODUCIBILITY.md](REPRODUCIBILITY.md). Full refits should run in a separate
copy because analysis scripts overwrite generated results. The historical
`run_review_simulation.R` uses the older optimization protocol; use the revised
optimization workflow in the reproduction guide for the current article.

The root `renv.lock` records the earlier simulation and movement environment.
It does not cover all genomic preprocessing dependencies. Their installation
instructions and recorded versions are provided separately in the genomic guide.

## Interpretation and limits

The genomic application concerns one cell line. Technical repeats and different
chromosomes do not supply biological replication. Three fitted allelic profiles
do not establish three biological populations, absolute copy-number states,
or a known true state count. The symmetry-aware F3 model has a lower BIC and
a higher confirmation score than M221; these comparisons are retained.
Nonrejection of the refinement constraints does not establish equivalence.

The supplementary elk results illustrate sensitivity to emission families and
temporal validation. They do not establish behavioral labels or a uniform
predictive advantage for mixture emissions.

## Compile and cite

With a LaTeX distribution installed, run `latexmk -pdf main.tex` from `article/`
and `latexmk -pdf supporting_information.tex` from `supplement/`.
Precompiled PDFs are included. Citation metadata are in [CITATION.cff](CITATION.cff).
Earlier public versions remain accessible through Git history.
