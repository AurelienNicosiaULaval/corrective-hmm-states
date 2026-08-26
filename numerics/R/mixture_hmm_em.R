# Expectation-maximization for univariate HMMs with Gaussian-mixture emissions.
# Source numerics/R/hmm_em.R before using this file.

component_mask <- function(component_counts, K, M) {
  component_counts <- as.integer(component_counts)
  if (length(component_counts) != K || anyNA(component_counts) ||
      any(component_counts < 1L) || any(component_counts > M)) {
    stop("component_counts must have K entries between 1 and M.", call. = FALSE)
  }
  outer(component_counts, seq_len(M), FUN = `>=`)
}

mixture_emission_probabilities <- function(y, weights, means, sds) {
  T <- length(y)
  K <- nrow(weights)
  M <- ncol(weights)
  component <- array(0, dim = c(T, K, M))
  for (state in seq_len(K)) {
    for (component_index in seq_len(M)) {
      component[, state, component_index] <- pmax(
        stats::dnorm(
          y,
          mean = means[state, component_index],
          sd = sds[state, component_index]
        ),
        1e-300
      )
    }
  }
  emission <- matrix(0, nrow = T, ncol = K)
  for (state in seq_len(K)) {
    component_matrix <- matrix(
      component[, state, ], nrow = T, ncol = M
    )
    emission[, state] <- component_matrix %*% weights[state, ]
  }
  list(emission = pmax(emission, 1e-300), component = component)
}

forward_backward_mixture <- function(
    y, initial, transition, weights, means, sds, lengths = NULL) {
  emission_result <- mixture_emission_probabilities(y, weights, means, sds)
  emission <- emission_result$emission
  K <- length(initial)
  M <- ncol(weights)
  posterior <- matrix(0, nrow = length(y), ncol = K)
  transition_sum <- matrix(0, nrow = K, ncol = K)
  log_likelihood <- 0

  for (indices in sequence_indices(length(y), lengths)) {
    result <- forward_backward_core(emission[indices, , drop = FALSE], initial, transition)
    posterior[indices, ] <- result$posterior
    transition_sum <- transition_sum + result$transition_sum
    log_likelihood <- log_likelihood + result$log_likelihood
  }

  component_posterior <- array(0, dim = c(length(y), K, M))
  for (state in seq_len(K)) {
    for (component_index in seq_len(M)) {
      component_posterior[, state, component_index] <-
        posterior[, state] * weights[state, component_index] *
        emission_result$component[, state, component_index] /
        emission[, state]
    }
  }

  list(
    log_likelihood = log_likelihood,
    posterior = posterior,
    transition_sum = transition_sum,
    component_posterior = component_posterior
  )
}

viterbi_single_mixture <- function(y, initial, transition, weights, means, sds) {
  T <- length(y)
  K <- length(initial)
  log_emission <- log(
    mixture_emission_probabilities(y, weights, means, sds)$emission
  )
  log_transition <- log(pmax(transition, 1e-300))
  delta <- matrix(0, nrow = T, ncol = K)
  back_pointer <- matrix(1L, nrow = T, ncol = K)
  delta[1L, ] <- log(pmax(initial, 1e-300)) + log_emission[1L, ]

  if (T > 1L) {
    for (time in 2:T) {
      for (destination in seq_len(K)) {
        candidates <- delta[time - 1L, ] + log_transition[, destination]
        back_pointer[time, destination] <- which.max(candidates)
        delta[time, destination] <- max(candidates) + log_emission[time, destination]
      }
    }
  }

  path <- integer(T)
  path[[T]] <- which.max(delta[T, ])
  if (T > 1L) {
    for (time in seq.int(T - 1L, 1L)) {
      path[[time]] <- back_pointer[time + 1L, path[[time + 1L]]]
    }
  }
  path
}

viterbi_mixture <- function(
    y, initial, transition, weights, means, sds, lengths = NULL) {
  path <- integer(length(y))
  for (indices in sequence_indices(length(y), lengths)) {
    path[indices] <- viterbi_single_mixture(
      y[indices], initial, transition, weights, means, sds
    )
  }
  path
}

initial_mixture_parameters <- function(y, K, M, component_counts, sigma_min) {
  y_sd <- max(stats::sd(y), sigma_min)
  state_centers <- as.numeric(
    stats::quantile(y, seq(0.15, 0.85, length.out = K), names = FALSE)
  ) + stats::rnorm(K, 0, 0.15 * y_sd)
  means <- matrix(rep(state_centers, M), nrow = K, ncol = M)
  sds <- matrix(max(2 * sigma_min, 0.45 * y_sd), nrow = K, ncol = M)
  weights <- matrix(0, nrow = K, ncol = M)

  for (state in seq_len(K)) {
    count <- component_counts[[state]]
    offsets <- if (count > 1L) seq(-0.35, 0.35, length.out = count) * y_sd else 0
    means[state, seq_len(count)] <- state_centers[[state]] + offsets +
      stats::rnorm(count, 0, 0.08 * y_sd)
    sds[state, seq_len(count)] <- sds[state, seq_len(count)] *
      stats::runif(count, 0.7, 1.3)
    weights[state, seq_len(count)] <- rdirichlet1(rep(2, count))
  }

  initial <- rdirichlet1(rep(1, K))
  transition <- t(vapply(
    seq_len(K),
    function(state) rdirichlet1(rep(1, K) + 6 * (seq_len(K) == state)),
    numeric(K)
  ))
  list(
    initial = initial,
    transition = transition,
    weights = weights,
    means = means,
    sds = sds
  )
}

sort_mixture_components <- function(weights, means, sds, component_counts) {
  for (state in seq_len(nrow(weights))) {
    count <- component_counts[[state]]
    ordering <- order(means[state, seq_len(count)])
    weights[state, seq_len(count)] <- weights[state, ordering]
    means[state, seq_len(count)] <- means[state, ordering]
    sds[state, seq_len(count)] <- sds[state, ordering]
  }
  list(weights = weights, means = means, sds = sds)
}

mixture_em_loop <- function(
    y, initial, transition, weights, means, sds, component_counts,
    max_iter, tolerance, sigma_min, lengths = NULL) {
  trace <- numeric()
  previous <- -Inf
  K <- length(initial)
  M <- ncol(weights)
  active <- component_mask(component_counts, K, M)
  sequences <- sequence_indices(length(y), lengths)
  starts <- vapply(sequences, function(index) index[[1L]], integer(1))

  for (iteration in seq_len(max_iter)) {
    fb <- forward_backward_mixture(
      y, initial, transition, weights, means, sds, lengths
    )
    log_likelihood <- fb$log_likelihood
    trace <- c(trace, log_likelihood)
    if (is.finite(previous) &&
        log_likelihood < previous - 1e-6 * (1 + abs(previous))) {
      stop(
        sprintf(
          "Mixture EM log-likelihood decreased at iteration %d: %.12f to %.12f.",
          iteration, previous, log_likelihood
        ),
        call. = FALSE
      )
    }
    converged <- is.finite(previous) &&
      (log_likelihood - previous) < tolerance * (1 + abs(previous))
    previous <- log_likelihood

    posterior <- fb$posterior
    initial <- normalize_probability(colMeans(posterior[starts, , drop = FALSE]) + 1e-10)
    transition <- normalize_rows(fb$transition_sum + 1e-10)

    component_totals <- apply(fb$component_posterior, c(2L, 3L), sum)
    weights <- ifelse(active, component_totals + 1e-12, 0)
    weights <- weights / rowSums(weights)
    safe_totals <- ifelse(active, component_totals + 1e-12, 1)

    updated_means <- matrix(0, nrow = K, ncol = M)
    updated_variances <- matrix(0, nrow = K, ncol = M)
    for (state in seq_len(K)) {
      for (component_index in seq_len(M)) {
        responsibility <- fb$component_posterior[, state, component_index]
        updated_means[state, component_index] <-
          sum(responsibility * y) / safe_totals[state, component_index]
      }
    }
    means[active] <- updated_means[active]
    for (state in seq_len(K)) {
      for (component_index in seq_len(M)) {
        responsibility <- fb$component_posterior[, state, component_index]
        updated_variances[state, component_index] <- sum(
          responsibility * (y - means[state, component_index])^2
        ) / safe_totals[state, component_index]
      }
    }
    sds[active] <- sqrt(pmax(updated_variances[active], sigma_min^2))

    sorted <- sort_mixture_components(weights, means, sds, component_counts)
    weights <- sorted$weights
    means <- sorted$means
    sds <- sorted$sds

    if (converged && iteration > 1L) {
      break
    }
  }

  fb <- forward_backward_mixture(
    y, initial, transition, weights, means, sds, lengths
  )
  trace <- c(trace, fb$log_likelihood)
  list(
    log_likelihood = fb$log_likelihood,
    initial = initial,
    transition = transition,
    weights = weights,
    means = means,
    sds = sds,
    posterior = fb$posterior,
    component_posterior = fb$component_posterior,
    n_iter = length(trace) - 1L,
    log_likelihood_trace = trace
  )
}

mixture_em_start <- function(
    y, K, M, component_counts, max_iter, tolerance, sigma_min,
    lengths = NULL) {
  initial_values <- initial_mixture_parameters(
    y, K, M, component_counts, sigma_min
  )
  mixture_em_loop(
    y,
    initial_values$initial,
    initial_values$transition,
    initial_values$weights,
    initial_values$means,
    initial_values$sds,
    component_counts,
    max_iter,
    tolerance,
    sigma_min,
    lengths
  )
}

mixture_n_parameters <- function(fit) {
  transition <- (fit$K - 1L) + fit$K * (fit$K - 1L)
  emissions <- sum((fit$component_counts - 1L) + 2L * fit$component_counts)
  as.integer(transition + emissions)
}

mixture_state_means <- function(fit) {
  rowSums(fit$weights * fit$means)
}

mixture_state_variances <- function(fit) {
  second <- rowSums(fit$weights * (fit$sds^2 + fit$means^2))
  second - mixture_state_means(fit)^2
}

fit_mixture_hmm <- function(
    y, K, M = 2L, seed = 1L, n_starts = 12L, max_iter = 250L,
    tolerance = 1e-6, sigma_min = SIGMA_MIN, screen_iter = NULL,
    refine_top = NULL, screen_tolerance = 1e-4,
    component_counts = rep(M, K), lengths = NULL) {
  K <- as.integer(K)
  M <- as.integer(M)
  component_counts <- as.integer(component_counts)
  component_mask(component_counts, K, M)
  set.seed(seed)

  if (!is.null(screen_iter) && !is.null(refine_top) && refine_top < n_starts) {
    screened <- lapply(seq_len(n_starts), function(unused) {
      mixture_em_start(
        y, K, M, component_counts, as.integer(screen_iter),
        screen_tolerance, sigma_min, lengths
      )
    })
    selected <- screened[order(
      vapply(screened, `[[`, numeric(1), "log_likelihood"),
      decreasing = TRUE
    )[seq_len(refine_top)]]
    candidates <- lapply(selected, function(result) {
      mixture_em_loop(
        y, result$initial, result$transition, result$weights,
        result$means, result$sds, component_counts, max_iter, tolerance,
        sigma_min, lengths
      )
    })
  } else {
    candidates <- lapply(seq_len(n_starts), function(unused) {
      mixture_em_start(
        y, K, M, component_counts, max_iter, tolerance, sigma_min, lengths
      )
    })
  }

  best <- candidates[[which.max(vapply(
    candidates, `[[`, numeric(1), "log_likelihood"
  ))]]
  state_means <- rowSums(best$weights * best$means)
  ordering <- order(state_means)
  fit <- list(
    model = "Mixture-emission HMM",
    K = K,
    M = M,
    component_counts = component_counts[ordering],
    log_likelihood = best$log_likelihood,
    initial = best$initial[ordering],
    transition = best$transition[ordering, ordering, drop = FALSE],
    weights = best$weights[ordering, , drop = FALSE],
    means = best$means[ordering, , drop = FALSE],
    sds = best$sds[ordering, , drop = FALSE],
    posterior = best$posterior[, ordering, drop = FALSE],
    component_posterior = best$component_posterior[, ordering, , drop = FALSE],
    n_iter = best$n_iter,
    log_likelihood_trace = best$log_likelihood_trace
  )
  fit$n_parameters <- mixture_n_parameters(fit)
  fit$viterbi <- viterbi_mixture(
    y, fit$initial, fit$transition, fit$weights, fit$means, fit$sds, lengths
  )
  class(fit) <- c("mixture_hmm_fit", "hmm_fit")
  fit
}

validate_mixture_fit <- function(fit, T) {
  active <- component_mask(fit$component_counts, fit$K, fit$M)
  stopifnot(
    inherits(fit, "mixture_hmm_fit"),
    identical(dim(fit$posterior), c(as.integer(T), fit$K)),
    identical(dim(fit$component_posterior), c(as.integer(T), fit$K, fit$M)),
    max(abs(rowSums(fit$transition) - 1)) < 1e-8,
    abs(sum(fit$initial) - 1) < 1e-8,
    max(abs(rowSums(fit$weights) - 1)) < 1e-8,
    all(abs(fit$weights[!active]) < 1e-12),
    all(fit$sds > 0),
    max(abs(rowSums(fit$posterior) - 1)) < 1e-8,
    all(is.finite(fit$posterior)),
    all(is.finite(fit$component_posterior)),
    all(is.finite(c(fit$log_likelihood, fit$means, fit$sds)))
  )
  differences <- diff(fit$log_likelihood_trace)
  tolerance <- -1e-6 * (1 + abs(head(fit$log_likelihood_trace, -1L)))
  stopifnot(all(differences >= tolerance))
  invisible(TRUE)
}
