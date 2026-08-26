# Expectation-maximization for a univariate Gaussian hidden Markov model.
#
# The implementation uses scaled forward-backward recursions, multiple random
# starts, a lower bound on emission standard deviations, and Viterbi decoding.
# No external HMM package is required.

SIGMA_MIN <- 0.05

normalize_probability <- function(x) {
  total <- sum(x)
  if (!is.finite(total) || total <= 0) {
    stop("Cannot normalize a non-positive probability vector.", call. = FALSE)
  }
  x / total
}

normalize_rows <- function(x) {
  totals <- rowSums(x)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop("Cannot normalize a matrix containing a non-positive row.", call. = FALSE)
  }
  x / totals
}

rdirichlet1 <- function(alpha) {
  normalize_probability(stats::rgamma(length(alpha), shape = alpha, rate = 1))
}

sequence_indices <- function(T, lengths = NULL) {
  if (is.null(lengths)) {
    return(list(seq_len(T)))
  }
  lengths <- as.integer(lengths)
  if (length(lengths) == 0L || anyNA(lengths) || any(lengths <= 0L)) {
    stop("lengths must contain positive sequence lengths.", call. = FALSE)
  }
  if (sum(lengths) != T) {
    stop("Sequence lengths must sum to the number of observations.", call. = FALSE)
  }
  ends <- cumsum(lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  Map(seq.int, starts, ends)
}

gaussian_emission_probabilities <- function(y, means, sds) {
  output <- vapply(
    seq_along(means),
    function(j) stats::dnorm(y, mean = means[[j]], sd = sds[[j]]),
    numeric(length(y))
  )
  output <- matrix(output, nrow = length(y), ncol = length(means))
  pmax(output, 1e-300)
}

forward_backward_core <- function(emission, initial, transition) {
  T <- nrow(emission)
  K <- ncol(emission)
  alpha <- matrix(0, nrow = T, ncol = K)
  scale <- numeric(T)

  alpha[1L, ] <- initial * emission[1L, ]
  scale[[1L]] <- sum(alpha[1L, ])
  alpha[1L, ] <- alpha[1L, ] / scale[[1L]]

  if (T > 1L) {
    for (time in 2:T) {
      alpha[time, ] <- as.numeric(alpha[time - 1L, ] %*% transition) *
        emission[time, ]
      scale[[time]] <- sum(alpha[time, ])
      alpha[time, ] <- alpha[time, ] / scale[[time]]
    }
  }

  beta <- matrix(0, nrow = T, ncol = K)
  beta[T, ] <- 1
  if (T > 1L) {
    for (time in seq.int(T - 1L, 1L)) {
      beta[time, ] <- as.numeric(
        transition %*% (emission[time + 1L, ] * beta[time + 1L, ])
      ) / scale[[time + 1L]]
    }
  }

  posterior <- alpha * beta
  posterior <- posterior / rowSums(posterior)
  transition_sum <- matrix(0, nrow = K, ncol = K)
  if (T > 1L) {
    for (time in seq_len(T - 1L)) {
      xi <- transition * outer(
        alpha[time, ],
        emission[time + 1L, ] * beta[time + 1L, ]
      )
      transition_sum <- transition_sum + xi / sum(xi)
    }
  }

  list(
    log_likelihood = sum(log(scale)),
    posterior = posterior,
    transition_sum = transition_sum
  )
}

forward_backward <- function(y, initial, transition, means, sds, lengths = NULL) {
  emission <- gaussian_emission_probabilities(y, means, sds)
  posterior <- matrix(0, nrow = length(y), ncol = length(initial))
  transition_sum <- matrix(0, nrow = length(initial), ncol = length(initial))
  log_likelihood <- 0

  for (indices in sequence_indices(length(y), lengths)) {
    result <- forward_backward_core(emission[indices, , drop = FALSE], initial, transition)
    posterior[indices, ] <- result$posterior
    transition_sum <- transition_sum + result$transition_sum
    log_likelihood <- log_likelihood + result$log_likelihood
  }

  list(
    log_likelihood = log_likelihood,
    posterior = posterior,
    transition_sum = transition_sum
  )
}

viterbi_single_gaussian <- function(y, initial, transition, means, sds) {
  T <- length(y)
  K <- length(initial)
  log_emission <- log(gaussian_emission_probabilities(y, means, sds))
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

viterbi_gaussian <- function(y, initial, transition, means, sds, lengths = NULL) {
  path <- integer(length(y))
  for (indices in sequence_indices(length(y), lengths)) {
    path[indices] <- viterbi_single_gaussian(
      y[indices], initial, transition, means, sds
    )
  }
  path
}

gaussian_em_loop <- function(
    y, initial, transition, means, sds, max_iter, tolerance, sigma_min,
    lengths = NULL) {
  trace <- numeric()
  previous <- -Inf
  indices <- sequence_indices(length(y), lengths)
  starts <- vapply(indices, function(index) index[[1L]], integer(1))

  for (iteration in seq_len(max_iter)) {
    fb <- forward_backward(y, initial, transition, means, sds, lengths)
    log_likelihood <- fb$log_likelihood
    trace <- c(trace, log_likelihood)

    if (is.finite(previous) &&
        log_likelihood < previous - 1e-6 * (1 + abs(previous))) {
      stop(
        sprintf(
          "EM log-likelihood decreased at iteration %d: %.12f to %.12f.",
          iteration, previous, log_likelihood
        ),
        call. = FALSE
      )
    }
    converged <- is.finite(previous) &&
      (log_likelihood - previous) < tolerance * (1 + abs(previous))
    previous <- log_likelihood

    posterior <- fb$posterior
    weights <- colSums(posterior) + 1e-12
    initial <- normalize_probability(colMeans(posterior[starts, , drop = FALSE]) + 1e-10)
    transition <- normalize_rows(fb$transition_sum + 1e-10)
    means <- colSums(posterior * y) / weights
    centered <- sweep(matrix(y, nrow = length(y), ncol = length(means)), 2L, means)
    variances <- colSums(posterior * centered^2) / weights
    sds <- sqrt(pmax(variances, sigma_min^2))

    if (converged && iteration > 1L) {
      break
    }
  }

  fb <- forward_backward(y, initial, transition, means, sds, lengths)
  trace <- c(trace, fb$log_likelihood)
  list(
    log_likelihood = fb$log_likelihood,
    initial = initial,
    transition = transition,
    means = means,
    sds = sds,
    posterior = fb$posterior,
    n_iter = length(trace) - 1L,
    log_likelihood_trace = trace
  )
}

gaussian_em_start <- function(
    y, K, max_iter, tolerance, sigma_min, lengths = NULL) {
  y_sd <- stats::sd(y)
  means <- as.numeric(stats::quantile(y, seq(0.1, 0.9, length.out = K), names = FALSE)) +
    stats::rnorm(K, 0, 0.25 * y_sd)
  sds <- rep(max(2 * sigma_min, y_sd / 2), K) * stats::runif(K, 0.6, 1.4)
  initial <- rdirichlet1(rep(1, K))
  transition <- t(vapply(
    seq_len(K),
    function(state) rdirichlet1(rep(1, K) + 6 * (seq_len(K) == state)),
    numeric(K)
  ))
  gaussian_em_loop(
    y, initial, transition, means, sds, max_iter, tolerance, sigma_min, lengths
  )
}

gaussian_n_parameters <- function(fit) {
  (fit$K - 1L) + fit$K * (fit$K - 1L) + 2L * fit$K
}

fit_bic <- function(fit, T) {
  -2 * fit$log_likelihood + fit$n_parameters * log(T)
}

fit_icl <- function(fit, T) {
  posterior <- pmax(fit$posterior, 1e-300)
  entropy <- -sum(fit$posterior * log(posterior))
  fit_bic(fit, T) + 2 * entropy
}

fit_gaussian_hmm <- function(
    y, K, seed, n_starts = 20L, max_iter = 300L, tolerance = 1e-6,
    sigma_min = SIGMA_MIN, screen_iter = NULL, refine_top = NULL,
    screen_tolerance = 1e-4, lengths = NULL) {
  stopifnot(length(K) == 1L, K >= 1L, n_starts >= 1L)
  set.seed(seed)

  if (!is.null(screen_iter) && !is.null(refine_top) && refine_top < n_starts) {
    screened <- lapply(seq_len(n_starts), function(unused) {
      gaussian_em_start(
        y, K, as.integer(screen_iter), screen_tolerance, sigma_min, lengths
      )
    })
    selected <- screened[order(
      vapply(screened, `[[`, numeric(1), "log_likelihood"),
      decreasing = TRUE
    )[seq_len(refine_top)]]
    candidates <- lapply(selected, function(result) {
      gaussian_em_loop(
        y, result$initial, result$transition, result$means, result$sds,
        max_iter, tolerance, sigma_min, lengths
      )
    })
  } else {
    candidates <- lapply(seq_len(n_starts), function(unused) {
      gaussian_em_start(y, K, max_iter, tolerance, sigma_min, lengths)
    })
  }

  best <- candidates[[which.max(vapply(
    candidates, `[[`, numeric(1), "log_likelihood"
  ))]]
  ordering <- order(best$means)
  fit <- list(
    model = "Gaussian HMM",
    K = as.integer(K),
    log_likelihood = best$log_likelihood,
    initial = best$initial[ordering],
    transition = best$transition[ordering, ordering, drop = FALSE],
    means = best$means[ordering],
    sds = best$sds[ordering],
    posterior = best$posterior[, ordering, drop = FALSE],
    n_iter = best$n_iter,
    log_likelihood_trace = best$log_likelihood_trace
  )
  fit$n_parameters <- gaussian_n_parameters(fit)
  fit$viterbi <- viterbi_gaussian(
    y, fit$initial, fit$transition, fit$means, fit$sds, lengths
  )
  class(fit) <- c("gaussian_hmm_fit", "hmm_fit")
  fit
}

validate_gaussian_fit <- function(fit, T) {
  stopifnot(
    inherits(fit, "gaussian_hmm_fit"),
    identical(dim(fit$posterior), c(as.integer(T), fit$K)),
    max(abs(rowSums(fit$transition) - 1)) < 1e-8,
    abs(sum(fit$initial) - 1) < 1e-8,
    all(fit$sds > 0),
    max(abs(rowSums(fit$posterior) - 1)) < 1e-8,
    all(is.finite(fit$posterior)),
    all(is.finite(c(fit$log_likelihood, fit$means, fit$sds)))
  )
  differences <- diff(fit$log_likelihood_trace)
  tolerance <- -1e-6 * (1 + abs(head(fit$log_likelihood_trace, -1L)))
  stopifnot(all(differences >= tolerance))
  invisible(TRUE)
}
