# True data-generating process used in the numerical study.

TRUE_TRANSITION <- matrix(
  c(0.93, 0.07,
    0.05, 0.95),
  nrow = 2,
  byrow = TRUE
)
TRUE_MIXTURE_WEIGHTS <- c(0.70, 0.30)
TRUE_REFINED_MEANS <- c(-1.20, 1.10, 3.80)
TRUE_REFINED_SDS <- c(0.55, 0.45, 0.65)
TRUE_REFINED_LABELS <- c("1a", "1b", "2")

stationary_distribution <- function(transition) {
  decomposition <- eigen(t(transition))
  index <- which.min(abs(decomposition$values - 1))
  stationary <- abs(Re(decomposition$vectors[, index]))
  stationary / sum(stationary)
}

simulate_hmm <- function(T, seed) {
  stopifnot(length(T) == 1L, T >= 1L)
  set.seed(seed)
  initial <- stationary_distribution(TRUE_TRANSITION)
  aggregate_state <- integer(T)
  refined_state <- integer(T)
  aggregate_state[[1L]] <- sample.int(2L, size = 1L, prob = initial)

  if (T > 1L) {
    for (time in 2:T) {
      aggregate_state[[time]] <- sample.int(
        2L,
        size = 1L,
        prob = TRUE_TRANSITION[aggregate_state[[time - 1L]], ]
      )
    }
  }

  for (time in seq_len(T)) {
    refined_state[[time]] <- if (aggregate_state[[time]] == 1L) {
      sample.int(2L, size = 1L, prob = TRUE_MIXTURE_WEIGHTS)
    } else {
      3L
    }
  }
  observations <- stats::rnorm(
    T,
    mean = TRUE_REFINED_MEANS[refined_state],
    sd = TRUE_REFINED_SDS[refined_state]
  )

  list(
    y = observations,
    aggregate_state = aggregate_state,
    refined_state = refined_state
  )
}
