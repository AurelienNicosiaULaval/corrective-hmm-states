test_that("stationary distribution is valid", {
  stationary <- stationary_distribution(TRUE_TRANSITION)
  expect_equal(as.numeric(stationary %*% TRUE_TRANSITION), stationary, tolerance = 1e-10)
  expect_equal(sum(stationary), 1, tolerance = 1e-12)
})

test_that("simulation returns coherent states", {
  simulation <- simulate_hmm(300, seed = 1)
  expect_length(simulation$y, 300)
  expect_length(simulation$aggregate_state, 300)
  expect_length(simulation$refined_state, 300)
  expect_true(all(simulation$aggregate_state %in% 1:2))
  expect_true(all(simulation$refined_state %in% 1:3))
  expect_identical(simulation$refined_state == 3L, simulation$aggregate_state == 2L)
  expect_true(all(is.finite(simulation$y)))
})

test_that("forward-backward posteriors and transition counts are normalized", {
  y <- simulate_hmm(200, seed = 2)$y
  result <- forward_backward(
    y,
    initial = c(0.5, 0.5),
    transition = matrix(c(0.9, 0.1, 0.1, 0.9), 2, byrow = TRUE),
    means = c(-1, 3),
    sds = c(1, 1)
  )
  expect_true(is.finite(result$log_likelihood))
  expect_equal(rowSums(result$posterior), rep(1, 200), tolerance = 1e-10)
  expect_equal(sum(result$transition_sum), 199, tolerance = 1e-6)
})

test_that("forward-backward handles independent sequences", {
  y <- simulate_hmm(220, seed = 22)$y
  result <- forward_backward(
    y,
    initial = c(0.5, 0.5),
    transition = matrix(c(0.9, 0.1, 0.1, 0.9), 2, byrow = TRUE),
    means = c(-1, 3),
    sds = c(1, 1),
    lengths = c(80, 60, 80)
  )
  expect_true(is.finite(result$log_likelihood))
  expect_equal(rowSums(result$posterior), rep(1, 220), tolerance = 1e-10)
  expect_equal(sum(result$transition_sum), 217, tolerance = 1e-6)
})

test_that("Gaussian EM fits are valid, ordered, and monotone", {
  y <- simulate_hmm(400, seed = 3)$y
  for (K in 2:3) {
    fit <- fit_gaussian_hmm(y, K, seed = 10 + K, n_starts = 4, max_iter = 150)
    expect_silent(validate_gaussian_fit(fit, length(y)))
    expect_equal(fit$K, K)
    expect_true(all(diff(fit$means) >= 0))
  }
})

test_that("BIC includes a positive parameter penalty", {
  y <- simulate_hmm(400, seed = 4)$y
  fit <- fit_gaussian_hmm(y, 2, seed = 5, n_starts = 3)
  expect_gt(fit_bic(fit, 400), -2 * fit$log_likelihood)
})

test_that("mixture-emission EM fit is valid, ordered, and monotone", {
  y <- simulate_hmm(260, seed = 6)$y
  fit <- fit_mixture_hmm(y, 2, M = 2, seed = 7, n_starts = 2, max_iter = 120)
  expect_silent(validate_mixture_fit(fit, length(y)))
  expect_equal(fit$K, 2)
  expect_equal(fit$M, 2)
  expect_true(all(diff(mixture_state_means(fit)) >= 0))
})

test_that("asymmetric components and multiple sequences are supported", {
  y <- simulate_hmm(260, seed = 26)$y
  fit <- fit_mixture_hmm(
    y,
    2,
    M = 2,
    component_counts = c(1, 2),
    lengths = c(130, 130),
    seed = 27,
    n_starts = 3,
    max_iter = 100
  )
  expect_silent(validate_mixture_fit(fit, length(y)))
  expect_equal(sort(fit$component_counts), c(1L, 2L))
  expect_equal(fit$n_parameters, 10)
  for (state in seq_len(fit$K)) {
    count <- fit$component_counts[[state]]
    if (count < fit$M) {
      expect_equal(fit$weights[state, seq.int(count + 1L, fit$M)], 0)
    }
  }
})

test_that("screened multi-start paths return valid fits", {
  y <- simulate_hmm(240, seed = 8)$y
  gaussian <- fit_gaussian_hmm(
    y, 3, seed = 9, n_starts = 5, max_iter = 60,
    screen_iter = 8, refine_top = 2
  )
  mixture <- fit_mixture_hmm(
    y, 2, M = 2, seed = 10, n_starts = 5, max_iter = 60,
    screen_iter = 8, refine_top = 2
  )
  expect_silent(validate_gaussian_fit(gaussian, length(y)))
  expect_silent(validate_mixture_fit(mixture, length(y)))
})

test_that("archived output files are present", {
  expected <- c(
    "numerics/output/review_gaussian_fit_summary.csv",
    "numerics/output/review_mixture_fit_summary.csv",
    "numerics/output/review_model_comparison.csv",
    "numerics/output/review_replications.csv",
    "numerics/output/review_summary_by_T.csv",
    "numerics/output/review_metadata.json",
    "article/figures/simulation_density.pdf",
    "article/figures/state_decoding_heatmap.pdf",
    "article/figures/bic_boxplot.pdf",
    "article/figures/bic_by_T.pdf",
    "article/figures/posterior_corrective_state.pdf",
    "article/figures/residual_acf.pdf",
    "article/tables/simulation_table.tex",
    "article/tables/replication_summary_table.tex",
    "empirical/output/elk_model_comparison.csv",
    "empirical/output/elk_decoding.csv",
    "empirical/figures/elk_model_diagnostics.pdf",
    "empirical/figures/elk_observation_diagnostics.pdf",
    "empirical/figures/elk_decoding_comparison.pdf",
    "empirical/figures/elk_state_mapping.pdf",
    "empirical/figures/elk_tracks_decoded.pdf",
    "article/tables/elk_model_comparison.tex",
    "supplement/tables/elk_mixture_parameters.tex",
    "supplement/tables/elk_predictive_diagnostics.tex"
  )
  paths <- file.path(ROOT, expected)
  expect_true(all(file.exists(paths)), info = paste(expected[!file.exists(paths)], collapse = "\n"))
})
