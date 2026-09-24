#!/usr/bin/env Rscript
# Four-state Gaussian comparator on the same chromosome-1 observations.
# Adds G4 without changing the previously fitted candidate models.
library(Rcpp)
library(data.table)
library(parallel)
library(digest)
source('numerics/R/hmm_em.R')
source('numerics/R/mixture_hmm_em.R')
source('numerics/R/revision_validation.R')
source('numerics/R/fast_em.R')
sourceCpp('numerics/R/forward_backward.cpp')
forward_backward_core <- forward_backward_compiled

input <- 'genomics/data/genomics/development_chr1.csv'
output <- 'genomics/results/genomic_pilot'
audit_output <- 'genomics/results/genomic_g4'
dir.create(audit_output, recursive = TRUE, showWarnings = FALSE)
observations <- fread(input)
y <- observations$beta_corrected
sequence_lengths <- as.integer(rle(observations$sequence)$lengths)
counts <- rep(1L, 4L)
stopifnot(all(observations$chromosome == 1L), all(is.finite(y)),
          sum(sequence_lengths) == length(y))
seed <- 20261314L
workers <- as.integer(Sys.getenv('GENOMIC_G4_WORKERS', '4'))
stopifnot(workers >= 1L)
reference <- readRDS(file.path(output, 'G3.rds'))
protocol <- list(seed = seed, n_random = 40L, n_refined = 8L,
                 screening_iterations = 100L, screening_tolerance = 1e-5,
                 refinement_iterations = 2500L, refinement_tolerance = 1e-8,
                 sigma_min = 0.01, sequence_lengths = sequence_lengths,
                 input_sha256 = digest(file = input, algo = 'sha256'),
                 reference_sha256 = digest(file = file.path(output, 'G3.rds'),
                                           algo = 'sha256'))
protocol_file <- file.path(audit_output, 'protocol.rds')
if (file.exists(protocol_file)) {
  stopifnot(identical(readRDS(protocol_file), protocol))
} else saveRDS(protocol, protocol_file)

fit_from <- function(initial_fit, max_iter, tolerance) {
  fast_em(y, initial_fit, counts, sequence_lengths,
          max_iter, tolerance, sigma_min = 0.01)
}
checkpoint <- function(label, initial_fit, max_iter, tolerance) {
  path <- file.path(audit_output, paste0(label, '.rds'))
  if (file.exists(path)) return(readRDS(path))
  fit <- fit_from(initial_fit, max_iter, tolerance)
  # Keep parameters and full likelihood histories for every start.
  fit$posterior <- NULL
  saveRDS(fit, path)
  cat(label, ': log likelihood =', format(fit$log_likelihood, digits = 12), '\n')
  fit
}

set.seed(seed)
initial_fits <- lapply(seq_len(40L), function(i) {
  initial_mixture_parameters(y, 4L, 1L, counts, 0.01)
})
screened <- mclapply(seq_along(initial_fits), function(i) {
  checkpoint(sprintf('screen_%02d', i), initial_fits[[i]], 100L, 1e-5)
}, mc.cores = workers)
if (any(vapply(screened, inherits, logical(1), 'try-error'))) stop('G4 screening failed')
screen_ll <- vapply(screened, `[[`, numeric(1), 'log_likelihood')
selected <- head(order(screen_ll, decreasing = TRUE), 8L)
refined <- mclapply(selected, function(i) {
  checkpoint(sprintf('refine_%02d', i), screened[[i]], 2500L, 1e-8)
}, mc.cores = workers)
if (any(vapply(refined, inherits, logical(1), 'try-error'))) stop('G4 refinement failed')

# Duplicating any G3 state yields an exactly equivalent feasible G4 start.
# Perturbing the duplicated means lets the two states subsequently separate.
nested <- list()
for (state in seq_len(3L)) {
  index <- c(seq_len(3L), state)
  split_weights <- rep(1, 4L)
  split_weights[c(state, 4L)] <- 0.5
  start <- list(initial = reference$initial[index] * split_weights,
                transition = sweep(reference$transition[index, index],
                                   2L, split_weights, '*'),
                weights = matrix(1, 4L, 1L),
                means = reference$means[index, , drop = FALSE],
                sds = reference$sds[index, , drop = FALSE])
  initial_ll <- forward_backward_mixture(y, start$initial, start$transition,
      start$weights, start$means, start$sds, sequence_lengths)$log_likelihood
  stopifnot(abs(initial_ll - reference$log_likelihood) < 1e-6)
  for (perturbation in c(0, 0.015, -0.015)) {
    candidate <- start
    candidate$means[c(state, 4L), 1L] <-
      candidate$means[c(state, 4L), 1L] + c(-perturbation, perturbation)
    label <- sprintf('nested_%d_%s', state,
                     c('zero', 'positive', 'negative')[match(perturbation, c(0, .015, -.015))])
    nested[[label]] <- checkpoint(label, candidate, 2500L, 1e-8)
  }
}
all_fits <- c(setNames(refined, sprintf('refine_%02d', selected)), nested)
candidate_ll <- vapply(all_fits, `[[`, numeric(1), 'log_likelihood')
best_label <- names(which.max(candidate_ll))
fit <- all_fits[[best_label]]
if (relative_gain(fit) > 1e-8) stop('Best G4 fit reached the iteration cap')
fit$K <- 4L
fit$M <- 1L
fit$component_counts <- counts
fit$n_parameters <- mixture_n_parameters(fit)
independent <- forward_backward_mixture(y, fit$initial, fit$transition,
    fit$weights, fit$means, fit$sds, sequence_lengths)
stopifnot(abs(independent$log_likelihood - fit$log_likelihood) < 1e-6)
fit$posterior <- independent$posterior
fit$viterbi <- viterbi_mixture(y, fit$initial, fit$transition,
    fit$weights, fit$means, fit$sds, sequence_lengths)
fit$BIC <- -2 * fit$log_likelihood + fit$n_parameters * log(length(y))
fit$converged <- TRUE
fit$selected_start <- best_label
fit$protocol <- protocol
fit$audit <- rbindlist(list(
  data.table(stage = 'screen', start = sprintf('screen_%02d', seq_along(screened)),
             log_likelihood = screen_ll),
  data.table(stage = 'refine', start = names(all_fits), log_likelihood = candidate_ll)))
stopifnot(fit$n_parameters == 23L,
          fit$log_likelihood >= reference$log_likelihood - 1e-6,
          readRDS(file.path(output, 'G5.rds'))$log_likelihood >= fit$log_likelihood - 1e-6)
saveRDS(fit, file.path(output, 'G4.rds'))
fwrite(fit$audit, file.path(audit_output, 'start_audit.csv'))

metrics <- fread(file.path(output, 'metrics.csv'))[model != 'G4']
new_row <- data.table(model = 'G4', states = 4L, components = 4L,
    parameters = fit$n_parameters, log_likelihood = fit$log_likelihood, BIC = fit$BIC,
    runs = sum(vapply(split(fit$viterbi, observations$sequence),
                      function(z) length(rle(z)$values), integer(1))),
    converged = fit$converged, relative_gain = relative_gain(fit),
    floor_hits = sum(fit$sds <= .01 * (1 + 1e-6)),
    minimum_transition = min(fit$transition))
fwrite(rbindlist(list(metrics, new_row)), file.path(output, 'metrics.csv'))
parameters <- fread(file.path(output, 'parameters.csv'))[model != 'G4']
new_parameters <- data.table(model = 'G4', state = seq_len(4L), component = 1L,
    weight = 1, mean = fit$means[, 1L], sd = fit$sds[, 1L],
    self_transition = diag(fit$transition), occupancy = colMeans(fit$posterior))
fwrite(rbindlist(list(parameters, new_parameters)), file.path(output, 'parameters.csv'))
observations[, state := fit$viterbi]
fwrite(observations, file.path(output, 'G4_decoded.csv'))
capture.output(sessionInfo(), file = file.path(audit_output, 'session_info.txt'))
print(new_row)
print(new_parameters)
