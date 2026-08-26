# Numerical summaries and observable diagnostics for fitted HMMs.

gaussian_fit_summary <- function(fit, T) {
  data.frame(
    K = fit$K,
    state = seq_len(fit$K),
    mean = fit$means,
    sd = fit$sds,
    posterior_occupation = colMeans(fit$posterior),
    logLik = fit$log_likelihood,
    BIC = fit_bic(fit, T),
    ICL = fit_icl(fit, T),
    check.names = FALSE
  )
}

mixture_fit_summary <- function(fit, T) {
  rows <- vector("list", fit$K * fit$M)
  index <- 1L
  occupations <- colMeans(fit$posterior)
  for (state in seq_len(fit$K)) {
    for (component in seq_len(fit$M)) {
      rows[[index]] <- data.frame(
        K = fit$K,
        M = fit$M,
        state = state,
        component = component,
        component_weight = fit$weights[state, component],
        mean = fit$means[state, component],
        sd = fit$sds[state, component],
        state_occupation = occupations[[state]],
        logLik = fit$log_likelihood,
        BIC = fit_bic(fit, T)
      )
      index <- index + 1L
    }
  }
  do.call(rbind, rows)
}

crosstab_true_fitted <- function(refined_state, fit, use = c("viterbi", "posterior")) {
  use <- match.arg(use)
  fitted <- if (use == "viterbi") fit$viterbi else max.col(fit$posterior)
  table(
    true_refined = factor(
      TRUE_REFINED_LABELS[refined_state],
      levels = TRUE_REFINED_LABELS
    ),
    fitted_state = factor(fitted, levels = seq_len(fit$K))
  )
}

one_step_residuals_gaussian <- function(y, fit) {
  prediction <- fit$initial
  residuals <- numeric(length(y))
  for (time in seq_along(y)) {
    prediction_mean <- sum(prediction * fit$means)
    prediction_second <- sum(prediction * (fit$sds^2 + fit$means^2))
    prediction_variance <- max(prediction_second - prediction_mean^2, 1e-10)
    residuals[[time]] <- (y[[time]] - prediction_mean) / sqrt(prediction_variance)
    density <- stats::dnorm(y[[time]], mean = fit$means, sd = fit$sds)
    filtered <- normalize_probability(prediction * pmax(density, 1e-300))
    prediction <- as.numeric(filtered %*% fit$transition)
  }
  residuals
}

residual_acf <- function(residuals, max_lag = 20L) {
  centered <- residuals - mean(residuals)
  denominator <- sum(centered^2)
  data.frame(
    lag = seq_len(max_lag),
    acf = vapply(seq_len(max_lag), function(lag) {
      if (denominator <= 0) {
        return(NA_real_)
      }
      sum(
        centered[seq_len(length(centered) - lag)] *
          centered[seq.int(lag + 1L, length(centered))]
      ) / denominator
    }, numeric(1))
  )
}
