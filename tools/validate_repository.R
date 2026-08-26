#!/usr/bin/env Rscript

# Validate the public reproducibility repository without refitting all models.

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
  "article/main.tex",
  "article/main.pdf",
  "article/references.bib",
  "article/figures/simulation_density.pdf",
  "article/figures/elk_observation_diagnostics.pdf",
  "article/figures/elk_state_mapping.pdf",
  "article/tables/elk_model_comparison.tex",
  "supplement/supporting_information.tex",
  "supplement/supporting_information.pdf",
  "empirical/DATA_SOURCE.md",
  "empirical/output/elk_analysis_metadata.json",
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
if (any(grepl("\\.py$", repository_files, ignore.case = TRUE))) {
  fail("Python source files remain in the R-only repository")
}
if (any(basename(repository_files) == "requirements.txt")) {
  fail("a Python requirements.txt file remains in the R-only repository")
}
if (any(tolower(basename(repository_files)) == "elk_data.csv")) {
  fail("elk_data.csv must not be distributed")
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
if (!identical(review_metadata$implementation, "R")) {
  fail("the archived simulation is not identified as an R implementation")
}
if (!identical(as.integer(review_metadata$n_replications_per_T), 200L)) {
  fail("archived simulation does not contain 200 replications per sample size")
}
if (!identical(as.integer(review_metadata$T_grid), c(500L, 1500L, 5000L))) {
  fail("unexpected simulation sample-size grid")
}

message("Required files: ", length(required_files), " found")
message("Implementation: R only")
message("Elk data: loaded from moveHMM; digest verified; no raw CSV distributed")
message("Archived simulation: 200 replications for T = 500, 1500, 5000")
message("Repository validation: passed")
