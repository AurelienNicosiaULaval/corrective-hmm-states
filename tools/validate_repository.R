#!/usr/bin/env Rscript

# Validate the public reproducibility repository without refitting all models.

library(jsonlite)
library(digest)

command_arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", command_arguments, value = TRUE)
if (length(script_argument) != 1L) {
  stop("Could not identify the validator path.", call. = FALSE)
}
script_path <- sub("^--file=", "", script_argument)
ROOT <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

EXPECTED_ELK_DIGEST <- "b18bb2b3a6f2b71bd7e95efcb59d9f5c907bd88a7664989e17ecea5ccbb1ea49"

required_files <- c(
  ".Rprofile",
  "DESCRIPTION",
  "renv.lock",
  "README.md",
  "CITATION.cff",
  "REPRODUCIBILITY.md",
  "manifest_sha256.json",
  "article/main.tex",
  "article/main.pdf",
  "article/references.bib",
  "article/figures/simulation_density.pdf",
  "article/figures/genomic_profiles.pdf",
  "article/figures/genomic_assays.pdf",
  "article/tables/genomic_comparison.tex",
  "genomics/README.md",
  "genomics/source_manifest.json",
  "genomics/scripts/test_genomic_results.R",
  "genomics/scripts/genomic_refinement_bootstrap.R",
  "genomics/results/genomic_pilot/metrics.csv",
  "genomics/results/refinement_bootstrap/summary.csv",
  "article/figures/elk_observation_diagnostics.pdf",
  "article/figures/elk_state_mapping.pdf",
  "article/tables/elk_model_comparison.tex",
  "supplement/supporting_information.tex",
  "supplement/supporting_information.pdf",
  "empirical/DATA_SOURCE.md",
  "empirical/output/elk_analysis_metadata.json",
  "empirical/output/elk_decoding.csv",
  "empirical/output/elk_model_comparison.csv",
  "numerics/R/hmm_em.R",
  "numerics/R/mixture_hmm_em.R",
  "numerics/R/plotting.R",
  "numerics/scripts/run_review_simulation.R",
  "numerics/scripts/run_elk_application.R",
  "numerics/output/review_metadata.json",
  "numerics/output/review_replications.csv",
  "numerics/tests/testthat.R"
)

fail <- function(message) {
  stop(paste0("Repository validation failed: ", message), call. = FALSE)
}

missing_files <- required_files[!file.exists(file.path(ROOT, required_files))]
if (length(missing_files) > 0L) {
  fail(paste("missing required files:", paste(missing_files, collapse = ", ")))
}

repository_files <- list.files(
  ROOT, recursive = TRUE, all.files = TRUE, full.names = TRUE,
  include.dirs = FALSE, no.. = TRUE
)
repository_files <- repository_files[!grepl("/(\\.git|renv/library|renv/staging)/", repository_files)]
# R/Rcpp implements the models; Python supports source downloads and extraction.
forbidden_elk_files <- c("elk_data.csv", "elk_daily_steps.csv")
if (any(tolower(basename(repository_files)) %in% forbidden_elk_files)) {
  fail("an elk coordinate or step-endpoint CSV is distributed")
}

empirical_csv_files <- list.files(
  file.path(ROOT, "empirical/output"),
  pattern = "\\.csv$",
  full.names = TRUE
)
coordinate_columns <- tolower(c(
  "Easting", "Northing", "x_from", "y_from", "x_to", "y_to"
))
offending_empirical_files <- vapply(empirical_csv_files, function(path) {
  header <- names(utils::read.csv(path, nrows = 0L, check.names = FALSE))
  any(tolower(header) %in% coordinate_columns)
}, logical(1))
if (any(offending_empirical_files)) {
  fail(paste(
    "coordinate endpoints occur in:",
    paste(basename(empirical_csv_files[offending_empirical_files]), collapse = ", ")
  ))
}

elk_metadata <- jsonlite::read_json(
  file.path(ROOT, "empirical/output/elk_analysis_metadata.json"),
  simplifyVector = TRUE
)
if (!identical(elk_metadata$implementation, "R")) {
  fail("the archived elk analysis is not identified as an R implementation")
}
if (!identical(elk_metadata$elk_data_source, "moveHMM::elk_data")) {
  fail("unexpected elk data source")
}
if (!identical(elk_metadata$elk_data_object_sha256, EXPECTED_ELK_DIGEST)) {
  fail("unexpected archived elk-data object digest")
}
if (!identical(elk_metadata$raw_data_csv_distributed, FALSE)) {
  fail("elk metadata does not state that raw coordinate CSVs are excluded")
}
expected_steps <- c("elk-115" = 193L, "elk-163" = 158L, "elk-287" = 163L, "elk-363" = 217L)
observed_steps <- unlist(elk_metadata$n_steps, use.names = TRUE)
if (!identical(as.integer(observed_steps[names(expected_steps)]), unname(expected_steps))) {
  fail("unexpected elk sample sizes")
}

elk_environment <- new.env(parent = emptyenv())
utils::data("elk_data", package = "moveHMM", envir = elk_environment)
if (!exists("elk_data", envir = elk_environment, inherits = FALSE)) {
  fail("moveHMM::elk_data could not be loaded")
}
elk_data <- get("elk_data", envir = elk_environment, inherits = FALSE)
observed_digest <- digest::digest(elk_data, algo = "sha256", serialize = TRUE)
if (!identical(observed_digest, EXPECTED_ELK_DIGEST)) {
  fail(paste("installed moveHMM::elk_data digest is", observed_digest))
}

review_metadata <- jsonlite::read_json(
  file.path(ROOT, "numerics/output/review_metadata.json"),
  simplifyVector = TRUE
)
if (!identical(review_metadata$implementation, "R with Rcpp EM verified against the R updates")) {
  fail("the archived simulation does not identify the revised R/Rcpp implementation")
}
if (!identical(as.integer(review_metadata$n_replications_per_T), 200L)) {
  fail("archived simulation does not contain 200 replications per sample size")
}
if (!identical(as.integer(review_metadata$T_grid), c(500L, 1500L, 5000L))) {
  fail("unexpected simulation sample-size grid")
}

manifest <- jsonlite::read_json(
  file.path(ROOT, "manifest_sha256.json"), simplifyVector = TRUE
)
if (anyDuplicated(manifest$path) || any(grepl("(^/|(^|/)\\.\\.(/|$))", manifest$path))) {
  fail("invalid or duplicated manifest paths")
}
for (i in seq_len(nrow(manifest))) {
  path <- file.path(ROOT, manifest$path[i])
  if (!file.exists(path) || file.info(path)$size != manifest$bytes[i] ||
      digest::digest(file = path, algo = "sha256") != manifest$sha256[i]) {
    fail(paste("manifest mismatch:", manifest$path[i]))
  }
}

expected_markers <- c(
  development_chr1 = 20697L, validation_chr12_17 = 52294L,
  confirmation_chr18_22 = 25413L, replication_chr1 = 20937L
)
for (name in names(expected_markers)) {
  path <- file.path(ROOT, "genomics/data/genomics", paste0(name, ".csv"))
  if (!file.exists(path) || nrow(utils::read.csv(path)) != expected_markers[[name]]) {
    fail(paste("unexpected genomic input size:", name))
  }
}
checkpoints <- list.files(
  file.path(ROOT, "numerics/output/optimization_audit"),
  pattern = "^fit_[0-9]+_[0-9]+\\.rds$"
)
bootstrap <- list.files(
  file.path(ROOT, "genomics/results/refinement_bootstrap"),
  pattern = "^replicate_[0-9]+\\.rds$"
)
if (length(checkpoints) != 600L || length(bootstrap) != 199L) {
  fail("incomplete corrected simulation or bootstrap checkpoints")
}

message("Required files: ", length(required_files), " found")
message("Distributed files: ", nrow(manifest), " SHA-256 digests verified")
message("Genomics: four processed inputs; 199 bootstrap checkpoints")
message("Implementation: R/Rcpp models; Python source-data utilities")
message("Elk data: loaded from moveHMM; digest verified; no coordinate endpoints distributed")
message("Archived simulation: 600 corrected checkpoints; 200 replications for T = 500, 1500, 5000")
message("Repository validation: passed")
