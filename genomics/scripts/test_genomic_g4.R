#!/usr/bin/env Rscript
# Independent checks of G4 and its retrospective comparisons.
library(Rcpp)
library(data.table)
library(testthat)
library(digest)
library(jsonlite)
source('numerics/R/hmm_em.R')
source('numerics/R/mixture_hmm_em.R')
source('numerics/R/revision_validation.R')
sourceCpp('numerics/R/forward_backward.cpp')
forward_backward_core <- forward_backward_compiled
sourceCpp('genomics/scripts/gaussian_mixture_hmm.cpp')

readfit <- function(name) readRDS(paste0('genomics/results/genomic_pilot/', name, '.rds'))
g4 <- readfit('G4')
development <- fread('genomics/data/genomics/development_chr1.csv')
lengths <- as.integer(rle(development$sequence)$lengths)
metrics <- fread('genomics/results/genomic_pilot/metrics.csv')

test_that('G4 likelihood, constraints and reported results agree independently', {
  expect_equal(g4$K, 4L)
  expect_equal(g4$component_counts, rep(1L, 4L))
  # Freely estimated initial law + transitions + Gaussian means and variances.
  expect_equal(g4$n_parameters, (4L - 1L) + 4L * (4L - 1L) + 2L * 4L)
  expect_true(all(is.finite(c(g4$means, g4$sds, g4$transition, g4$initial))))
  expect_true(all(g4$transition >= 0))
  expect_true(all(g4$initial >= 0))
  expect_equal(rowSums(g4$transition), rep(1, 4), tolerance = 1e-12)
  expect_equal(sum(g4$initial), 1, tolerance = 1e-12)
  expect_true(all(g4$sds > 0.01))
  expect_lte(relative_gain(g4), 1e-8)
  expect_true(all(diff(g4$log_likelihood_trace) >= -1e-6))
  evaluated <- forward_backward_mixture(development$beta_corrected, g4$initial,
      g4$transition, g4$weights, g4$means, g4$sds, lengths)
  expect_equal(evaluated$log_likelihood, g4$log_likelihood, tolerance = 1e-11)
  expect_equal(evaluated$posterior, g4$posterior, tolerance = 1e-10)
  expect_equal(g4$BIC, -2 * evaluated$log_likelihood + 23 * log(nrow(development)),
               tolerance = 1e-11)
  expect_gte(g4$log_likelihood, readfit('G3')$log_likelihood)
  expect_gte(readfit('G5')$log_likelihood, g4$log_likelihood)
  expect_equal(metrics[model == 'G4', log_likelihood], g4$log_likelihood)
  expect_equal(metrics[model == 'G4', BIC], g4$BIC)
  parameters <- fread('genomics/results/genomic_pilot/parameters.csv')[model == 'G4']
  expect_equal(parameters$occupancy, colMeans(evaluated$posterior), tolerance = 1e-10)
  decoded <- fread('genomics/results/genomic_pilot/G4_decoded.csv')
  expect_equal(decoded$state, g4$viterbi)
  expect_equal(metrics[model == 'G4', runs],
    sum(vapply(split(decoded$state, decoded$sequence), function(z) length(rle(z)$values), integer(1))))
})

test_that('screened and nested starts remain auditable and selection is by likelihood', {
  audit <- fread('genomics/results/genomic_g4/start_audit.csv')
  expect_equal(sum(audit$stage == 'screen'), 40L)
  expect_equal(sum(audit$stage == 'refine'), 17L)
  expect_equal(sum(startsWith(audit$start, 'nested_')), 9L)
  expect_equal(g4$log_likelihood, max(audit[stage == 'refine', log_likelihood]))
  expect_equal(g4$log_likelihood, audit[start == g4$selected_start, log_likelihood])
  expect_identical(g4$protocol$input_sha256,
                   digest(file = 'genomics/data/genomics/development_chr1.csv', algo = 'sha256'))
  expect_identical(g4$protocol$reference_sha256,
                   digest(file = 'genomics/results/genomic_pilot/G3.rds', algo = 'sha256'))
})

test_that('all G4 scores are independently recomputed on the stated fixed inputs', {
  sources <- c(replication = 'replication_chr1', validation = 'validation_chr12_17',
               confirmation = 'confirmation_chr18_22')
  saved <- fread('genomics/results/genomic_validation/scores.csv')
  by_chr <- fread('genomics/results/genomic_validation/chromosome_scores.csv')
  for (dataset_name in names(sources)) {
    d <- fread(paste0('genomics/data/genomics/', sources[[dataset_name]], '.csv'))
    starts <- c(TRUE, diff(d$sequence) != 0)
    emissions <- sapply(seq_len(4), function(k) dnorm(d$beta_corrected, g4$means[k, 1], g4$sds[k, 1]))
    evaluated <- fb_sequences_cpp(pmax(emissions, 1e-300), g4$initial, g4$transition, starts)
    summary <- saved[model == 'G4' & dataset == dataset_name]
    expect_equal(nrow(summary), 1L)
    expect_equal(summary$n, nrow(d))
    expect_equal(summary$log_score, evaluated$log_likelihood, tolerance = 1e-10)
    predictions <- fread(paste0('genomics/results/genomic_validation/', dataset_name, '_G4.csv'))
    expect_equal(predictions$log_score, log(evaluated$scales), tolerance = 1e-10)
    expect_equal(predictions$state,
                 viterbi_sequences_cpp(log(pmax(emissions, 1e-300)), g4$initial, g4$transition, starts))
    expect_equal(sum(by_chr[model == 'G4' & dataset == dataset_name, log_score]),
                 evaluated$log_likelihood, tolerance = 1e-10)
  }
  # Retain the comparison favourable to G4, including its adverse later block.
  expect_gt(saved[model == 'G4' & dataset == 'validation', log_score],
            saved[model == 'M221' & dataset == 'validation', log_score])
  expect_lt(saved[model == 'G4' & dataset == 'confirmation', log_score],
            saved[model == 'M221' & dataset == 'confirmation', log_score])
})

test_that('the seven original model objects retain their pre-evaluation hashes', {
  frozen <- read_json('genomics/audit/confirmation_freeze.json')$files
  models <- names(frozen)[grepl('\\.rds$', names(frozen))]
  expect_length(models, 7L)
  expect_false(any(grepl('G4', models)))
  for (path in models) expect_identical(digest(file = file.path('genomics', path), algo = 'sha256'), frozen[[path]])
})
cat('G4 and retrospective-comparison checks completed.\n')
