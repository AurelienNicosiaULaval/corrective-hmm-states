# Targeted second-revision analyses. Source after the original EM functions.

relative_gain <- function(fit) {
  tail(diff(fit$log_likelihood_trace), 1) /
    (1 + abs(tail(fit$log_likelihood_trace, 2)[1]))
}

refine_candidate <- function(y, fit, counts, floor, iterations, tolerance) {
  mixture_em_loop(y, fit$initial, fit$transition, fit$weights, fit$means,
                 fit$sds, counts, iterations, tolerance, floor)
}

# The expanded transition matrix is the exact refinement of Proposition 2.5.
expand_mixture_fit <- function(fit) {
  state <- rep(seq_len(fit$K), fit$component_counts)
  component <- unlist(lapply(fit$component_counts, seq_len))
  index <- cbind(state, component)
  w <- fit$weights[index]
  transition <- fit$transition[state, state, drop = FALSE] *
    matrix(w, length(w), length(w), byrow = TRUE)
  list(initial = fit$initial[state] * w, transition = transition,
       weights = matrix(1, length(w), 1), means = matrix(fit$means[index], ncol=1),
       sds = matrix(fit$sds[index], ncol=1))
}

fit_audited <- function(y, counts, floor, seed, nested = NULL,
                        n_starts = 200L, refine_top = 20L) {
  set.seed(seed)
  K <- length(counts); M <- max(counts)
  screened <- lapply(seq_len(n_starts), function(i)
    mixture_em_start(y, K, M, counts, 70L, 1e-4, floor))
  ranking <- order(vapply(screened, `[[`, numeric(1), "log_likelihood"), decreasing=TRUE)
  selected <- ranking[seq_len(refine_top)]
  candidates <- lapply(screened[selected], function(f)
    refine_candidate(y, f, counts, floor, 2000L, 1e-8))
  if (!is.null(nested)) {
    candidates <- c(candidates, list(refine_candidate(y, expand_mixture_fit(nested),
                                                       counts, floor, 2000L, 1e-8)))
    selected <- c(selected, 0L)
  }
  # Continue any unfinished refined run, keeping its entire likelihood trace.
  candidates <- lapply(candidates, function(f) {
    for (extra in seq_len(4L)) {
      if (relative_gain(f) <= 1e-8) break
      updated <- refine_candidate(y, f, counts, floor, 2000L, 1e-8)
      updated$log_likelihood_trace <- c(f$log_likelihood_trace,
                                        updated$log_likelihood_trace[-1])
      updated$n_iter <- f$n_iter + updated$n_iter
      f <- updated
    }
    f
  })
  audit_rows <- function(items, stage, index) do.call(rbind, lapply(seq_along(items), function(i) {
    f <- items[[i]]
    data.frame(stage=stage, start=index[i], log_likelihood=f$log_likelihood,
               iterations=f$n_iter, relative_gain=relative_gain(f),
               converged=relative_gain(f)<=1e-8,
               minimum_transition=min(f$transition),
               maximum_self_transition=max(diag(f$transition)),
               floor_hits=sum(f$sds <= floor*(1+1e-6)))
  }))
  best_index <- which.max(vapply(candidates, `[[`, numeric(1), "log_likelihood"))
  fit <- candidates[[best_index]]
  ordering <- order(rowSums(fit$weights*fit$means))
  fit$K <- as.integer(K); fit$M <- as.integer(M)
  fit$component_counts <- as.integer(counts[ordering])
  fit$initial <- fit$initial[ordering]
  fit$transition <- fit$transition[ordering,ordering,drop=FALSE]
  fit$weights <- fit$weights[ordering,,drop=FALSE]
  fit$means <- fit$means[ordering,,drop=FALSE]
  fit$sds <- fit$sds[ordering,,drop=FALSE]
  fit$posterior <- fit$posterior[,ordering,drop=FALSE]
  fit$component_posterior <- fit$component_posterior[,ordering,,drop=FALSE]
  fit$n_parameters <- mixture_n_parameters(fit)
  fit$viterbi <- viterbi_mixture(y,fit$initial,fit$transition,fit$weights,fit$means,fit$sds)
  fit$audit <- rbind(audit_rows(screened,"screen",seq_len(n_starts)),
                     audit_rows(candidates,"refine",selected))
  fit$best_start <- selected[best_index]
  fit$seed <- seed; fit$sigma_min <- floor
  fit$converged <- relative_gain(fit) <= 1e-8
  class(fit) <- c("mixture_hmm_fit","hmm_fit")
  validate_mixture_fit(fit,length(y))
  fit
}

# Parameters stay fixed at the training estimate. Filtering sees only past y.
sequential_predictions <- function(y, fit) {
  pred <- fit$initial
  out <- matrix(NA_real_,length(y),3,
                dimnames=list(NULL,c("log_score","pit","negative_mass")))
  for (t in seq_along(y)) {
    dens <- cdf <- negative <- numeric(fit$K)
    for (k in seq_len(fit$K)) {
      j <- seq_len(fit$component_counts[k])
      dens[k] <- sum(fit$weights[k,j]*stats::dnorm(y[t],fit$means[k,j],fit$sds[k,j]))
      cdf[k] <- sum(fit$weights[k,j]*stats::pnorm(y[t],fit$means[k,j],fit$sds[k,j]))
      negative[k] <- sum(fit$weights[k,j]*stats::pnorm(0,fit$means[k,j],fit$sds[k,j]))
    }
    out[t,] <- c(log(max(sum(pred*dens),1e-300)),sum(pred*cdf),sum(pred*negative))
    pred <- as.numeric(normalize_probability(pred*pmax(dens,1e-300)) %*% fit$transition)
  }
  as.data.frame(out)
}

density_decomposition <- function(fit, grid) {
  occupation <- colMeans(fit$posterior)
  output <- list()
  for (k in seq_len(fit$K)) for (j in seq_len(fit$component_counts[k])) {
    output[[length(output)+1L]] <- data.frame(y=grid,state=k,component=j,
      density=occupation[k]*fit$weights[k,j]*stats::dnorm(grid,fit$means[k,j],fit$sds[k,j]))
  }
  do.call(rbind,output)
}

count_modes <- function(fit) {
  # Include all component tails; modes refer to the real-line fitted density.
  grid <- seq(min(fit$means-5*fit$sds),max(fit$means+5*fit$sds),length.out=50001L)
  parts <- density_decomposition(fit,grid)
  density <- rowsum(parts$density,rep(seq_along(grid),sum(fit$component_counts)))
  modes <- which(diff(sign(diff(as.vector(density)))) == -2L)+1L
  list(count=length(modes),locations=grid[modes])
}

canonical_parameters <- function(fit, counts) {
  ordering <- order(fit$component_counts, decreasing=TRUE)
  stopifnot(identical(as.integer(fit$component_counts[ordering]),as.integer(counts)))
  list(initial=fit$initial[ordering],transition=fit$transition[ordering,ordering,drop=FALSE],
       weights=fit$weights[ordering,,drop=FALSE],means=fit$means[ordering,,drop=FALSE],
       sds=fit$sds[ordering,,drop=FALSE])
}

adopt_candidate <- function(fit, candidate, counts, y, stage) {
  row <- data.frame(stage=stage,start=0L,log_likelihood=candidate$log_likelihood,
    iterations=candidate$n_iter,relative_gain=relative_gain(candidate),
    converged=relative_gain(candidate)<=1e-8,minimum_transition=min(candidate$transition),
    maximum_self_transition=max(diag(candidate$transition)),
    floor_hits=sum(candidate$sds<=fit$sigma_min*(1+1e-6)))
  fit$audit <- rbind(fit$audit,row)
  if(candidate$log_likelihood > fit$log_likelihood+1e-7) {
    ordering <- order(rowSums(candidate$weights*candidate$means))
    fit$log_likelihood <- candidate$log_likelihood
    fit$log_likelihood_trace <- candidate$log_likelihood_trace
    fit$n_iter <- candidate$n_iter
    fit$initial <- candidate$initial[ordering]
    fit$transition <- candidate$transition[ordering,ordering,drop=FALSE]
    for(name in c("weights","means","sds")) fit[[name]] <- candidate[[name]][ordering,,drop=FALSE]
    fit$posterior <- candidate$posterior[,ordering,drop=FALSE]
    fit$component_posterior <- candidate$component_posterior[,ordering,,drop=FALSE]
    fit$component_counts <- as.integer(counts[ordering])
    fit$viterbi <- viterbi_mixture(y,fit$initial,fit$transition,fit$weights,fit$means,fit$sds)
    fit$converged <- relative_gain(fit)<=1e-8
    fit$best_start <- stage
  }
  validate_mixture_fit(fit,length(y))
  fit
}
