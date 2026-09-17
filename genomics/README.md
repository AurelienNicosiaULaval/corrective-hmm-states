# Genomic application reproduction

Run all commands from the parent reproducibility-package directory. Analyses use R 4.5.0 and a C++ compiler supported by Rcpp. The recorded package inventory and processing session are in `audit`. The original study is [Chiang et al. (2009), Nature Methods 6, 99–103](https://doi.org/10.1038/nmeth.1276); the four selected arrays are from [GEO GSE13372](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13372). See the [input table in the repository README](../README.md#where-are-the-data) and [source manifest](source_manifest.json) for the supplied files, original download URLs and SHA-256 hashes.

## Reproduce results from supplied inputs

The package contains the full exported signals, hg18 coordinates, filtered marker tables, fitted objects, diagnostic tables and source hashes. The HCC1143 arrays measure one tumor cell line with a matched normal cell line; the repeated arrays are technical repeats. No ground-truth allelic-state labels are available.

Before genomic commands, run `export RENV_CONFIG_AUTOLOADER_ENABLED=false` in the terminal to disable the legacy project autoloader.

Core requirements are Rcpp, data.table, ggplot2, patchwork, clue and testthat. To reproduce preprocessing and the domain comparators, also install aroma.affymetrix, aroma.light, affxparser, PSCBS, DNAcopy and sfit. `scripts/install_dependencies.R` installs into the isolated `genomics/library` directory. It targets Bioconductor 3.21 and the recorded sfit commit; the package inventory documents versions used for the reported results. Set `R_LIBS_USER` to the absolute path of this isolated library when needed. The package-level historical renv lock belongs to the retained simulation and elk work; it is not a lock for the new genomic dependencies.

```
Rscript genomics/scripts/test_genomic_results.R
Rscript genomics/scripts/test_symmetric_genomic_hmm.R
Rscript genomics/scripts/evaluate_genomic_pilot.R
Rscript genomics/scripts/genomic_diagnostics.R
```

The evaluation script uses fixed fitted parameters and resets filtering at recorded sequence boundaries. The file named `confirmation` contains only chromosomes 18–22 and represents a later held-out block in an exploratory analysis. The seven fitted objects and the original protocol, preserved in `audit/protocol_before_confirmation.md`, match the hashes in `audit/confirmation_freeze.json`. The current `GENOMICS_PROTOCOL.md` is an edited reader guide; `audit/protocol_provenance.json` records that distinction. The distributed evaluation script is a later version and does not match its earlier recorded hash. This is not a preregistered confirmatory study. Keep the saved models fixed when reproducing the reported scores.

## Complete refit

Use a separate working copy, since the following commands overwrite generated fits and result tables.

```
Rscript genomics/scripts/genomic_pilot.R
Rscript genomics/scripts/genomic_symmetry_benchmark.R
Rscript genomics/scripts/genomic_sensitivity.R
Rscript genomics/scripts/genomic_pscbs_comparison.R
Rscript genomics/scripts/evaluate_genomic_pilot.R
Rscript genomics/scripts/genomic_diagnostics.R
Rscript genomics/scripts/test_genomic_results.R
```

The first command refits all five unrestricted Gaussian/mixture candidates from fixed random seeds, including exact nesting and persistent variance-regime starts. F2/F3 are symmetry-aware comparators. Model choice is always by likelihood within a specified family. Numerical precision can vary across compilers; the independent refit on the reported platform reproduced the selected likelihoods and decoded paths.

## Regenerate from original microarrays

This route requires several gigabytes of local space and the larger preprocessing dependencies. It is unnecessary for reproducing the reported HMM estimates from the supplied processed inputs.

```
python3 genomics/scripts/fetch_sources.py
Rscript genomics/scripts/process_snp6.R
Rscript genomics/scripts/prepare_genomic_development.R
Rscript genomics/scripts/prepare_genomic_validation.R
Rscript genomics/scripts/prepare_genomic_confirmation.R
python3 genomics/scripts/extract_sequencing_table.py
```

The Python extraction command requires pypdf. Source downloads are checked against SHA-256 hashes. Annotation hashes concern their uncompressed representation. TumorBoost is applied after selecting normal candidate heterozygotes; corrected values outside [0,1] are retained. The processing script exports hg18 coordinates for analysis and independent sequencing comparison. The hg19 fragment-length annotation used during standard preprocessing does not change the hg18 coordinates in the analysis tables.

The full raw preprocessing was performed for the analysis. The delivered raw-regeneration instructions fix the coordinate export that was initially performed as a separate step. The supplied processed-input path, numerical refit, predictive evaluation, table extraction and tests were verified separately; the raw download-and-preprocess sequence was not repeated from a second empty machine.

## Interpretation and limitations

The application illustrates a statistical modeling problem in one cell line. The fitted states describe regional allelic profiles; they do not identify cell populations, absolute copy-number states, or a known true order. F3 has lower BIC than M221 and a higher confirmation score. PSCBS uses total intensity as well as fractions. Published sequencing segments provide an independent assay of total copy ratio, not ground-truth allelic-state labels. Broad assignments are repeatable, but short runs are sensitive to marker thinning and residual magnitude dependence remains.

## Refinement constraints

Run `Rscript genomics/scripts/genomic_refinement_bootstrap.R` for the
model-conditional parametric bootstrap comparing M221 with G5. The script
resumes saved replicates and retains all fitted parameters and likelihood
traces; see `results/refinement_bootstrap/protocol.txt`.
