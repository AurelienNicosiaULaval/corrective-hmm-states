#!/usr/bin/env Rscript

# Reproduce the numerical study, derived tables, and ggplot2 figures.
# Run from the reproducibility-package root:
#   Rscript numerics/scripts/run_review_simulation.R

required_packages <- c("dplyr", "ggplot2", "jsonlite")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop(
    "Install the required R packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

if (!file.exists("numerics/R/hmm_em.R")) {
  stop("Run this script from the reproducibility-package root.", call. = FALSE)
}

source("numerics/R/hmm_em.R")
source("numerics/R/mixture_hmm_em.R")
source("numerics/R/simulate.R")
source("numerics/R/diagnostics.R")
source("numerics/R/plotting.R")

environment_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = as.character(default))
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  parsed
}

SEED <- 20260610L
T_MAIN <- 1500L
T_GRID <- c(500L, 1500L, 5000L)
N_REPS <- environment_integer("REVIEW_N_REPS", 200L)
N_STARTS_MAIN <- environment_integer("REVIEW_N_STARTS_MAIN", 10L)
N_STARTS_REPS <- environment_integer("REVIEW_N_STARTS_REPS", 10L)
N_STARTS_MIX_MAIN <- environment_integer("REVIEW_N_STARTS_MIX_MAIN", 10L)
N_STARTS_MIX_REPS <- environment_integer("REVIEW_N_STARTS_MIX_REPS", 10L)
MAX_ITER_REPS <- environment_integer("REVIEW_MAX_ITER_REPS", 100L)
SCREEN_ITER_REPS <- environment_integer("REVIEW_SCREEN_ITER_REPS", 25L)
REFINE_TOP_REPS <- environment_integer("REVIEW_REFINE_TOP_REPS", 3L)
SCREEN_ITER_MIX_REPS <- environment_integer("REVIEW_SCREEN_ITER_MIX_REPS", 30L)
REFINE_TOP_MIX_REPS <- environment_integer("REVIEW_REFINE_TOP_MIX_REPS", 3L)
N_WORKERS <- environment_integer(
  "REVIEW_N_WORKERS",
  max(1L, min(8L, parallel::detectCores(logical = FALSE) - 1L))
)

OUTPUT_DIRECTORY <- "numerics/output"
FIGURE_DIRECTORY <- "article/figures"
ARTICLE_TABLE_DIRECTORY <- "article/tables"
SUPPLEMENT_TABLE_DIRECTORY <- "supplement/tables"
PARTIAL_RESULTS_PATH <- file.path(OUTPUT_DIRECTORY, "review_replications.partial.csv")
PARTIAL_METADATA_PATH <- file.path(OUTPUT_DIRECTORY, "review_replications.partial.json")

write_csv <- function(data, path, row.names = FALSE) {
  utils::write.csv(data, path, row.names = row.names, na = "")
}

fit_replication <- function(task) {
  T <- task$T
  replication <- task$rep
  seed <- SEED + 1000000L + 10000L * T + replication
  y <- simulate_hmm(T, seed)$y
  gaussian2 <- fit_gaussian_hmm(
    y, 2, seed = seed + 2L, n_starts = N_STARTS_REPS,
    max_iter = MAX_ITER_REPS, screen_iter = SCREEN_ITER_REPS,
    refine_top = REFINE_TOP_REPS
  )
  gaussian3 <- fit_gaussian_hmm(
    y, 3, seed = seed + 3L, n_starts = N_STARTS_REPS,
    max_iter = MAX_ITER_REPS, screen_iter = SCREEN_ITER_REPS,
    refine_top = REFINE_TOP_REPS
  )
  mixture2 <- fit_mixture_hmm(
    y, 2, M = 2, seed = seed + 20L, n_starts = N_STARTS_MIX_REPS,
    max_iter = MAX_ITER_REPS, screen_iter = SCREEN_ITER_MIX_REPS,
    refine_top = REFINE_TOP_MIX_REPS
  )
  validate_gaussian_fit(gaussian2, T)
  validate_gaussian_fit(gaussian3, T)
  validate_mixture_fit(mixture2, T)

  data.frame(
    T = T,
    rep = replication,
    seed = seed,
    logLik_gauss_K2 = gaussian2$log_likelihood,
    logLik_gauss_K3 = gaussian3$log_likelihood,
    logLik_mix_K2_M2 = mixture2$log_likelihood,
    BIC_gauss_K2 = fit_bic(gaussian2, T),
    BIC_gauss_K3 = fit_bic(gaussian3, T),
    BIC_mix_K2_M2 = fit_bic(mixture2, T),
    delta_BIC_gauss_K3_minus_K2 = fit_bic(gaussian3, T) - fit_bic(gaussian2, T),
    delta_BIC_mix_K2_minus_gauss_K3 = fit_bic(mixture2, T) - fit_bic(gaussian3, T),
    gaussian_BIC_favours_K3 = fit_bic(gaussian3, T) < fit_bic(gaussian2, T),
    enriched_BIC_favours_K2_mix = fit_bic(mixture2, T) < fit_bic(gaussian3, T),
    mix_state1_mean = mixture_state_means(mixture2)[[1L]],
    mix_state2_mean = mixture_state_means(mixture2)[[2L]],
    mix_state1_occ = colMeans(mixture2$posterior)[[1L]],
    mix_state2_occ = colMeans(mixture2$posterior)[[2L]]
  )
}

partial_configuration <- function() {
  list(
    implementation = "R",
    seed = SEED,
    T_grid = T_GRID,
    n_replications_per_T = N_REPS,
    n_starts_replications = N_STARTS_REPS,
    n_starts_mixture_replications = N_STARTS_MIX_REPS,
    max_iter_replications = MAX_ITER_REPS,
    screen_iter_replications = SCREEN_ITER_REPS,
    refine_top_replications = REFINE_TOP_REPS,
    screen_iter_mixture_replications = SCREEN_ITER_MIX_REPS,
    refine_top_mixture_replications = REFINE_TOP_MIX_REPS
  )
}

configuration_json <- function(configuration) {
  jsonlite::toJSON(configuration, auto_unbox = TRUE, digits = 16, pretty = TRUE)
}

load_partial_results <- function() {
  configuration <- partial_configuration()
  if (!file.exists(PARTIAL_RESULTS_PATH) || !file.exists(PARTIAL_METADATA_PATH)) {
    writeLines(configuration_json(configuration), PARTIAL_METADATA_PATH)
    return(data.frame())
  }
  existing <- paste(readLines(PARTIAL_METADATA_PATH, warn = FALSE), collapse = "\n")
  if (!identical(existing, configuration_json(configuration))) {
    message("Ignoring partial results because their configuration differs.")
    unlink(PARTIAL_RESULTS_PATH)
    writeLines(configuration_json(configuration), PARTIAL_METADATA_PATH)
    return(data.frame())
  }
  output <- utils::read.csv(PARTIAL_RESULTS_PATH, stringsAsFactors = FALSE)
  output <- output[output$T %in% T_GRID & output$rep >= 1L & output$rep <= N_REPS, ]
  output <- output[!duplicated(output[c("T", "rep")], fromLast = TRUE), ]
  output[order(output$T, output$rep), ]
}

write_partial_results <- function(results) {
  writeLines(configuration_json(partial_configuration()), PARTIAL_METADATA_PATH)
  if (nrow(results) == 0L) {
    return(invisible(NULL))
  }
  results <- results[!duplicated(results[c("T", "rep")], fromLast = TRUE), ]
  results <- results[order(results$T, results$rep), ]
  temporary <- paste0(PARTIAL_RESULTS_PATH, ".tmp")
  write_csv(results, temporary)
  if (!file.rename(temporary, PARTIAL_RESULTS_PATH)) {
    stop("Could not update the partial replication file.", call. = FALSE)
  }
  invisible(NULL)
}

run_replications <- function(start_time) {
  task_grid <- expand.grid(
    T = sort(T_GRID, decreasing = TRUE),
    rep = seq_len(N_REPS),
    KEEP.OUT.ATTRS = FALSE
  )
  tasks <- split(task_grid, seq_len(nrow(task_grid)))
  results <- load_partial_results()
  completed <- if (nrow(results) == 0L) character() else paste(results$T, results$rep)
  remaining <- tasks[!vapply(
    tasks, function(task) paste(task$T, task$rep) %in% completed, logical(1)
  )]
  message(
    sprintf(
      "Running %d remaining replications out of %d with %d worker(s).",
      length(remaining), length(tasks), N_WORKERS
    )
  )

  if (length(remaining) == 0L) {
    return(results)
  }
  batches <- split(remaining, ceiling(seq_along(remaining) / 25L))
  for (batch in batches) {
    batch_results <- if (.Platform$OS.type == "unix" && N_WORKERS > 1L) {
      parallel::mclapply(batch, fit_replication, mc.cores = N_WORKERS, mc.preschedule = FALSE)
    } else {
      lapply(batch, fit_replication)
    }
    results <- rbind(results, do.call(rbind, batch_results))
    write_partial_results(results)
    counts <- table(factor(results$T, levels = T_GRID))
    progress <- paste(sprintf("T=%d: %d/%d", T_GRID, counts, N_REPS), collapse = ", ")
    message(sprintf(
      "%d/%d replications complete [%s] (%.0f s elapsed).",
      nrow(results), length(tasks), progress,
      as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    ))
  }
  results
}

model_rows <- function(gaussian2, gaussian3, mixture2, T) {
  data.frame(
    model = c("Gaussian HMM", "Gaussian HMM", "Mixture-emission HMM"),
    emission = c(
      "one Gaussian per state",
      "one Gaussian per state",
      "two Gaussians per state"
    ),
    K = c(2L, 3L, 2L),
    logLik = c(
      gaussian2$log_likelihood,
      gaussian3$log_likelihood,
      mixture2$log_likelihood
    ),
    BIC = c(
      fit_bic(gaussian2, T),
      fit_bic(gaussian3, T),
      fit_bic(mixture2, T)
    ),
    diagnostic = c(
      "merges the heterogeneous regime",
      "splits the heterogeneous regime",
      "absorbs the split in the emission"
    )
  )
}

write_article_model_table <- function(rows) {
  specification <- c(
    "Gaussian, one Gaussian/state",
    "Gaussian, one Gaussian/state",
    "Mixture, two Gaussians/state"
  )
  interpretation <- c("merged regime", "split regime", "emission split")
  lines <- c(
    "% Generated by numerics/scripts/run_review_simulation.R.",
    "\\begin{table}[t]",
    "\\centering",
    "\\small",
    "\\caption{Main numerical diagnostic, $T=1500$. Lower BIC is better.}",
    "\\label{tab:sim}",
    "\\begin{tabular}{@{}lrrrp{0.24\\linewidth}@{}}",
    "\\toprule",
    "Specification & $K$ & logLik & BIC & interpretation\\\\",
    "\\midrule"
  )
  for (index in seq_len(nrow(rows))) {
    lines <- c(lines, sprintf(
      "%s & %d & %.1f & %.1f & %s\\\\",
      specification[[index]], rows$K[[index]], rows$logLik[[index]],
      rows$BIC[[index]], interpretation[[index]]
    ))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, file.path(ARTICLE_TABLE_DIRECTORY, "simulation_table.tex"))
}

write_replication_tables <- function(summary) {
  lines <- c(
    "% Generated by numerics/scripts/run_review_simulation.R.",
    "\\begin{table}[t]",
    "\\centering",
    "\\small",
    paste0(
      "\\caption{BIC selection frequencies over fixed-seed replications. Lower BIC is better; ",
      "$\\Delta_{\\rm G}=\\mathrm{BIC}_{G,3}-\\mathrm{BIC}_{G,2}$ and ",
      "$\\Delta_{\\rm M}=\\mathrm{BIC}_{M,2}-\\mathrm{BIC}_{G,3}$.}"
    ),
    "\\label{tab:replications}",
    "\\begin{tabular}{rrrrr}",
    "\\toprule",
    "$T$ & reps & mean $\\Delta_{\\rm G}$ & Pr$(\\Delta_{\\rm G}<0)$ & Pr$(\\Delta_{\\rm M}<0)$\\\\",
    "\\midrule"
  )
  for (index in seq_len(nrow(summary))) {
    lines <- c(lines, sprintf(
      "%d & %d & %.1f & %.3f & %.3f\\\\",
      summary$T[[index]], summary$n_replications[[index]],
      summary$mean_delta_BIC_gauss_K3_minus_K2[[index]],
      summary$prop_gaussian_BIC_favours_K3[[index]],
      summary$prop_enriched_BIC_favours_K2_mix[[index]]
    ))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, file.path(ARTICLE_TABLE_DIRECTORY, "replication_summary_table.tex"))
  writeLines(lines, file.path(SUPPLEMENT_TABLE_DIRECTORY, "review_summary_table.tex"))
}

main <- function() {
  start_time <- Sys.time()
  invisible(lapply(
    c(OUTPUT_DIRECTORY, FIGURE_DIRECTORY, ARTICLE_TABLE_DIRECTORY, SUPPLEMENT_TABLE_DIRECTORY),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))

  message("Fitting the reference-series models.")
  simulation <- simulate_hmm(T_MAIN, SEED)
  y <- simulation$y
  gaussian2 <- fit_gaussian_hmm(
    y, 2, seed = SEED + 2L, n_starts = N_STARTS_MAIN
  )
  gaussian3 <- fit_gaussian_hmm(
    y, 3, seed = SEED + 3L, n_starts = N_STARTS_MAIN
  )
  mixture2 <- fit_mixture_hmm(
    y, 2, M = 2, seed = SEED + 20L, n_starts = N_STARTS_MIX_MAIN
  )
  validate_gaussian_fit(gaussian2, T_MAIN)
  validate_gaussian_fit(gaussian3, T_MAIN)
  validate_mixture_fit(mixture2, T_MAIN)

  gaussian_summary <- rbind(
    gaussian_fit_summary(gaussian2, T_MAIN),
    gaussian_fit_summary(gaussian3, T_MAIN)
  )
  write_csv(gaussian_summary, file.path(OUTPUT_DIRECTORY, "review_gaussian_fit_summary.csv"))
  write_csv(
    mixture_fit_summary(mixture2, T_MAIN),
    file.path(OUTPUT_DIRECTORY, "review_mixture_fit_summary.csv")
  )

  table2 <- crosstab_true_fitted(simulation$refined_state, gaussian2)
  table3 <- crosstab_true_fitted(simulation$refined_state, gaussian3)
  write_csv(as.data.frame.matrix(table2), file.path(OUTPUT_DIRECTORY, "crosstab_fit2.csv"), row.names = TRUE)
  write_csv(as.data.frame.matrix(table3), file.path(OUTPUT_DIRECTORY, "crosstab_fit3.csv"), row.names = TRUE)

  plot_simulation_density(y, gaussian2, gaussian3, mixture2, FIGURE_DIRECTORY)
  plot_decoding_heatmap(table2, table3, FIGURE_DIRECTORY)
  corrective_state <- which.min(abs(gaussian3$means - TRUE_REFINED_MEANS[[2L]]))
  plot_corrective_posterior(y, gaussian3, corrective_state, FIGURE_DIRECTORY)
  acf2 <- residual_acf(one_step_residuals_gaussian(y, gaussian2))
  acf3 <- residual_acf(one_step_residuals_gaussian(y, gaussian3))
  write_csv(acf2, file.path(OUTPUT_DIRECTORY, "residual_acf_fit2.csv"))
  write_csv(acf3, file.path(OUTPUT_DIRECTORY, "residual_acf_fit3.csv"))
  plot_residual_acf(acf2, acf3, FIGURE_DIRECTORY)
  message("Reference-series models and figures are complete.")

  replications <- run_replications(start_time)
  replications <- replications[order(replications$T, replications$rep), ]
  rownames(replications) <- NULL
  write_csv(replications, file.path(OUTPUT_DIRECTORY, "review_replications.csv"))
  summary <- dplyr::summarise(
    dplyr::group_by(replications, T),
    n_replications = dplyr::n(),
    mean_delta_BIC_gauss_K3_minus_K2 = mean(delta_BIC_gauss_K3_minus_K2),
    prop_gaussian_BIC_favours_K3 = mean(gaussian_BIC_favours_K3),
    mean_delta_BIC_mix_K2_minus_gauss_K3 = mean(delta_BIC_mix_K2_minus_gauss_K3),
    prop_enriched_BIC_favours_K2_mix = mean(enriched_BIC_favours_K2_mix),
    .groups = "drop"
  )
  write_csv(summary, file.path(OUTPUT_DIRECTORY, "review_summary_by_T.csv"))
  plot_bic_boxplot(
    replications$delta_BIC_gauss_K3_minus_K2[replications$T == T_MAIN],
    FIGURE_DIRECTORY
  )
  plot_bic_by_sample_size(summary, FIGURE_DIRECTORY)

  comparison <- model_rows(gaussian2, gaussian3, mixture2, T_MAIN)
  write_csv(comparison, file.path(OUTPUT_DIRECTORY, "review_model_comparison.csv"))
  write_article_model_table(comparison)
  write_replication_tables(summary)

  metadata <- list(
    implementation = "R",
    r_version = R.version.string,
    seed = SEED,
    T_main = T_MAIN,
    T_grid = T_GRID,
    n_replications_per_T = N_REPS,
    n_starts_main = N_STARTS_MAIN,
    n_starts_replications = N_STARTS_REPS,
    n_starts_mixture_main = N_STARTS_MIX_MAIN,
    n_starts_mixture_replications = N_STARTS_MIX_REPS,
    n_workers = N_WORKERS,
    max_iter_replications = MAX_ITER_REPS,
    screen_iter_replications = SCREEN_ITER_REPS,
    refine_top_replications = REFINE_TOP_REPS,
    screen_iter_mixture_replications = SCREEN_ITER_MIX_REPS,
    refine_top_mixture_replications = REFINE_TOP_MIX_REPS,
    true_transition = split(TRUE_TRANSITION, row(TRUE_TRANSITION)),
    true_stationary = stationary_distribution(TRUE_TRANSITION),
    true_mixture_weights_regime1 = TRUE_MIXTURE_WEIGHTS,
    true_refined_means = TRUE_REFINED_MEANS,
    true_refined_sds = TRUE_REFINED_SDS,
    main_model_comparison = unname(split(comparison, seq_len(nrow(comparison)))),
    selection_summary = unname(split(summary, seq_len(nrow(summary)))),
    runtime_seconds = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  )
  jsonlite::write_json(
    metadata,
    file.path(OUTPUT_DIRECTORY, "review_metadata.json"),
    auto_unbox = TRUE,
    digits = 16,
    pretty = TRUE
  )
  unlink(c(PARTIAL_RESULTS_PATH, PARTIAL_METADATA_PATH))
  print(comparison)
  print(summary)
}

main()
