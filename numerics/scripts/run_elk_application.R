#!/usr/bin/env Rscript

# Empirical diagnostic application to Canadian elk movement data.
#
# The data are read directly from moveHMM::elk_data. No CSV containing the
# original coordinates or derived coordinate endpoints is written or
# redistributed. Run from the reproducibility-package root:
#   Rscript numerics/scripts/run_elk_application.R

required_packages <- c(
  "digest", "dplyr", "ggplot2", "jsonlite", "moveHMM", "patchwork", "scales", "tidyr"
)
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
source("numerics/R/plotting.R")

environment_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value) || value < 1L) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  value
}

SEED <- 20260826L
N_STARTS <- environment_integer("ELK_N_STARTS", 200L)
SCREEN_ITER <- environment_integer("ELK_SCREEN_ITER", 70L)
REFINE_TOP <- environment_integer("ELK_REFINE_TOP", 10L)
MAX_ITER <- environment_integer("ELK_MAX_ITER", 500L)
MINIMUM_SD <- 0.05
EXPECTED_DATA_DIGEST <- "b18bb2b3a6f2b71bd7e95efcb59d9f5c907bd88a7664989e17ecea5ccbb1ea49"

OUTPUT_DIRECTORY <- "empirical/output"
FIGURE_DIRECTORY <- "empirical/figures"
ARTICLE_FIGURE_DIRECTORY <- "article/figures"
ARTICLE_TABLE_DIRECTORY <- "article/tables"
SUPPLEMENT_TABLE_DIRECTORY <- "supplement/tables"

MODEL_LABELS <- c(
  gaussian_k2 = "Gaussian K=2",
  gaussian_k3 = "Gaussian K=3",
  mixture_k2 = "Mixture K=2"
)
MODEL_COLORS <- c(
  gaussian_k2 = "#C76E2E",
  gaussian_k3 = "#2F6B9A",
  mixture_k2 = "#2D8A66"
)
STATE_COLORS <- c("1" = "#4C78A8", "2" = "#F2A541", "3" = "#C44E52")

load_elk_data <- function() {
  environment <- new.env(parent = emptyenv())
  utils::data("elk_data", package = "moveHMM", envir = environment)
  if (!exists("elk_data", envir = environment, inherits = FALSE)) {
    stop("moveHMM::elk_data could not be loaded.", call. = FALSE)
  }
  elk_data <- get("elk_data", envir = environment, inherits = FALSE)
  expected_names <- c("ID", "Easting", "Northing", "dist_water")
  if (!identical(names(elk_data), expected_names) || nrow(elk_data) != 735L) {
    stop("moveHMM::elk_data does not have the expected structure.", call. = FALSE)
  }
  observed_digest <- digest::digest(elk_data, algo = "sha256", serialize = TRUE)
  if (!identical(observed_digest, EXPECTED_DATA_DIGEST)) {
    stop(
      "The installed moveHMM::elk_data object does not match the archived analysis input. ",
      "Expected SHA-256 ", EXPECTED_DATA_DIGEST, ", observed ", observed_digest, ".",
      call. = FALSE
    )
  }
  elk_data
}

prepare_steps <- function(raw_data) {
  tracks <- split(raw_data, raw_data$ID, drop = TRUE)
  output <- lapply(tracks, function(track) {
    track <- track[seq_len(nrow(track)), , drop = FALSE]
    displacement_x <- diff(as.numeric(track$Easting))
    displacement_y <- diff(as.numeric(track$Northing))
    step_km <- sqrt(displacement_x^2 + displacement_y^2) / 1000
    data.frame(
      ID = as.character(track$ID[-nrow(track)]),
      step_index = seq_along(step_km),
      x_from = as.numeric(track$Easting[-nrow(track)]) / 1000,
      y_from = as.numeric(track$Northing[-nrow(track)]) / 1000,
      x_to = as.numeric(track$Easting[-1L]) / 1000,
      y_to = as.numeric(track$Northing[-1L]) / 1000,
      step_km = step_km,
      log1p_step_km = log1p(step_km),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, unname(output))
}

fit_elk_models <- function(y, animal_index) {
  common <- list(
    n_starts = N_STARTS,
    max_iter = MAX_ITER,
    screen_iter = SCREEN_ITER,
    refine_top = REFINE_TOP,
    sigma_min = MINIMUM_SD
  )
  gaussian2 <- do.call(
    fit_gaussian_hmm,
    c(list(y = y, K = 2L, seed = SEED + 1000L * animal_index + 2L), common)
  )
  gaussian3 <- do.call(
    fit_gaussian_hmm,
    c(list(y = y, K = 3L, seed = SEED + 1000L * animal_index + 3L), common)
  )
  mixture2 <- do.call(
    fit_mixture_hmm,
    c(
      list(
        y = y, K = 2L, M = 2L, component_counts = c(2L, 2L),
        seed = SEED + 1000L * animal_index + 20L
      ),
      common
    )
  )
  validate_gaussian_fit(gaussian2, length(y))
  validate_gaussian_fit(gaussian3, length(y))
  validate_mixture_fit(mixture2, length(y))
  list(gaussian_k2 = gaussian2, gaussian_k3 = gaussian3, mixture_k2 = mixture2)
}

pit_gaussian <- function(y, fit) {
  prediction <- fit$initial
  output <- numeric(length(y))
  for (time in seq_along(y)) {
    output[[time]] <- sum(
      prediction * stats::pnorm(y[[time]], mean = fit$means, sd = fit$sds)
    )
    emission <- stats::dnorm(y[[time]], mean = fit$means, sd = fit$sds)
    filtered <- normalize_probability(prediction * pmax(emission, 1e-300))
    prediction <- as.numeric(filtered %*% fit$transition)
  }
  pmin(pmax(output, 1e-10), 1 - 1e-10)
}

pit_mixture <- function(y, fit) {
  prediction <- fit$initial
  output <- numeric(length(y))
  for (time in seq_along(y)) {
    state_cdf <- numeric(fit$K)
    state_density <- numeric(fit$K)
    for (state in seq_len(fit$K)) {
      components <- seq_len(fit$component_counts[[state]])
      state_cdf[[state]] <- sum(
        fit$weights[state, components] * stats::pnorm(
          y[[time]],
          mean = fit$means[state, components],
          sd = fit$sds[state, components]
        )
      )
      state_density[[state]] <- sum(
        fit$weights[state, components] * stats::dnorm(
          y[[time]],
          mean = fit$means[state, components],
          sd = fit$sds[state, components]
        )
      )
    }
    output[[time]] <- sum(prediction * state_cdf)
    filtered <- normalize_probability(prediction * pmax(state_density, 1e-300))
    prediction <- as.numeric(filtered %*% fit$transition)
  }
  pmin(pmax(output, 1e-10), 1 - 1e-10)
}

pit_metrics <- function(values) {
  sorted <- sort(values)
  n <- length(sorted)
  ks <- max(max(seq_len(n) / n - sorted), max(sorted - (seq_len(n) - 1L) / n))
  cvm <- 1 / (12 * n) + sum((sorted - (2 * seq_len(n) - 1) / (2 * n))^2)
  coverage <- mean(sorted >= 0.05 & sorted <= 0.95)
  c(PIT_KS_distance = ks, PIT_CvM = cvm, PIT_90pct_coverage = coverage)
}

model_row <- function(animal, y, model_key, fit) {
  pit <- if (model_key == "mixture_k2") pit_mixture(y, fit) else pit_gaussian(y, fit)
  metrics <- pit_metrics(pit)
  data.frame(
    ID = animal,
    n_steps = length(y),
    model_key = model_key,
    model = unname(MODEL_LABELS[[model_key]]),
    logLik = fit$log_likelihood,
    n_parameters = fit$n_parameters,
    BIC = fit_bic(fit, length(y)),
    PIT_KS_distance = metrics[["PIT_KS_distance"]],
    PIT_CvM = metrics[["PIT_CvM"]],
    PIT_90pct_coverage = metrics[["PIT_90pct_coverage"]]
  )
}

parameter_rows <- function(animal, fits) {
  gaussian_rows <- list()
  row_index <- 1L
  for (model_key in c("gaussian_k2", "gaussian_k3")) {
    fit <- fits[[model_key]]
    for (state in seq_len(fit$K)) {
      gaussian_rows[[row_index]] <- data.frame(
        ID = animal,
        model = unname(MODEL_LABELS[[model_key]]),
        state = state,
        mean_log1p_step = fit$means[[state]],
        median_step_km = expm1(fit$means[[state]]),
        sd_log1p_step = fit$sds[[state]],
        posterior_occupation = mean(fit$posterior[, state]),
        self_transition = fit$transition[state, state]
      )
      row_index <- row_index + 1L
    }
  }

  fit <- fits$mixture_k2
  mixture_rows <- list()
  row_index <- 1L
  for (state in seq_len(fit$K)) {
    for (component in seq_len(fit$component_counts[[state]])) {
      mixture_rows[[row_index]] <- data.frame(
        ID = animal,
        model = unname(MODEL_LABELS[["mixture_k2"]]),
        state = state,
        state_label = paste("aggregate", state),
        component = component,
        component_weight = fit$weights[state, component],
        mean_log1p_step = fit$means[state, component],
        median_step_km = expm1(fit$means[state, component]),
        sd_log1p_step = fit$sds[state, component],
        posterior_state_occupation = mean(fit$posterior[, state]),
        self_transition = fit$transition[state, state]
      )
      row_index <- row_index + 1L
    }
  }
  list(
    gaussian = do.call(rbind, gaussian_rows),
    mixture = do.call(rbind, mixture_rows)
  )
}

mapping_rows <- function(animal, gaussian3, mixture2) {
  counts <- table(
    factor(gaussian3$viterbi, levels = 1:3),
    factor(mixture2$viterbi, levels = 1:2)
  )
  proportions <- prop.table(counts, margin = 1)
  output <- expand.grid(
    gaussian_k3_state = 1:3,
    mixture_k2_state = 1:2,
    KEEP.OUT.ATTRS = FALSE
  )
  output$ID <- animal
  output$count <- as.vector(counts)
  output$row_proportion <- as.vector(proportions)
  output[c("ID", "gaussian_k3_state", "mixture_k2_state", "count", "row_proportion")]
}

transition_rows <- function(animal, fits) {
  rows <- list()
  row_index <- 1L
  for (model_key in names(fits)) {
    fit <- fits[[model_key]]
    for (origin in seq_len(fit$K)) {
      for (destination in seq_len(fit$K)) {
        rows[[row_index]] <- data.frame(
          ID = animal,
          model = unname(MODEL_LABELS[[model_key]]),
          from_state = origin,
          to_state = destination,
          probability = fit$transition[origin, destination]
        )
        row_index <- row_index + 1L
      }
    }
  }
  do.call(rbind, rows)
}

write_model_table <- function(model_table) {
  pivot <- tidyr::pivot_wider(
    model_table[c("ID", "n_steps", "model_key", "BIC")],
    names_from = model_key,
    values_from = BIC
  )
  lines <- c(
    "% Generated by numerics/scripts/run_elk_application.R.",
    "\\begin{table}[t]",
    "\\centering",
    "\\small",
    paste0(
      "\\caption{Empirical comparison for four Canadian elk. Lower BIC is better. ",
      "The enriched model has two aggregate states and two Gaussian components ",
      "in each state-dependent distribution.}"
    ),
    "\\label{tab:elk-comparison}",
    "\\begin{tabular}{lrrrrr}",
    "\\toprule",
    "Animal & steps & Gaussian $K=2$ & Gaussian $K=3$ & Mixture $K=2$ & $\\Delta$BIC\\\\",
    "\\midrule"
  )
  for (index in seq_len(nrow(pivot))) {
    delta <- pivot$mixture_k2[[index]] - pivot$gaussian_k3[[index]]
    lines <- c(lines, sprintf(
      "%s & %d & %.1f & %.1f & %.1f & %+.1f\\\\",
      pivot$ID[[index]], pivot$n_steps[[index]], pivot$gaussian_k2[[index]],
      pivot$gaussian_k3[[index]], pivot$mixture_k2[[index]], delta
    ))
  }
  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{minipage}{0.92\\linewidth}\\footnotesize",
    paste0(
      "\\textit{Note:} $\\Delta$BIC is mixture $K=2$ minus Gaussian $K=3$; ",
      "negative values favour the two-state enriched-emission model."
    ),
    "\\end{minipage}",
    "\\end{table}",
    ""
  )
  writeLines(lines, file.path(ARTICLE_TABLE_DIRECTORY, "elk_model_comparison.tex"))
}

write_parameter_table <- function(mixture_table) {
  lines <- c(
    "% Generated by numerics/scripts/run_elk_application.R.",
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    paste0(
      "\\caption{State-dependent parameters for the enriched two-state elk models. ",
      "Component medians are back-transformed to kilometres.}"
    ),
    "\\label{tab:elk-mixture-parameters}",
    "\\begin{tabular}{llrrrr@{\\hspace{2em}}r}",
    "\\toprule",
    "Animal & state & comp. & weight & median & occupation & self-transition\\\\",
    "\\midrule"
  )
  for (index in seq_len(nrow(mixture_table))) {
    row <- mixture_table[index, ]
    lines <- c(lines, sprintf(
      "%s & %s & %d & %.2f & %.2f km & %.2f & %.2f\\\\",
      row$ID, row$state_label, row$component, row$component_weight,
      row$median_step_km, row$posterior_state_occupation, row$self_transition
    ))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}", "")
  writeLines(lines, file.path(SUPPLEMENT_TABLE_DIRECTORY, "elk_mixture_parameters.tex"))
}

write_predictive_table <- function(model_table) {
  lines <- c(
    "% Generated by numerics/scripts/run_elk_application.R.",
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    paste0(
      "\\caption{One-step-ahead predictive calibration for the elk models. Smaller ",
      "Kolmogorov--Smirnov (KS) distance from uniformity is better; nominal interval ",
      "coverage is 0.90.}"
    ),
    "\\label{tab:elk-predictive}",
    "\\begin{tabular}{llrr}",
    "\\toprule",
    "Animal & model & PIT KS distance & 90\\% coverage\\\\",
    "\\midrule"
  )
  for (index in seq_len(nrow(model_table))) {
    row <- model_table[index, ]
    lines <- c(lines, sprintf(
      "%s & %s & %.3f & %.3f\\\\",
      row$ID, row$model, row$PIT_KS_distance, row$PIT_90pct_coverage
    ))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}", "")
  writeLines(lines, file.path(SUPPLEMENT_TABLE_DIRECTORY, "elk_predictive_diagnostics.tex"))
}

save_elk_plot <- function(plot, name, width, height) {
  save_ggplot(plot, FIGURE_DIRECTORY, name, width, height)
  ggplot2::ggsave(
    file.path(ARTICLE_FIGURE_DIRECTORY, paste0(name, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )
  invisible(plot)
}

plot_model_diagnostics <- function(model_table) {
  baseline <- model_table[model_table$model_key == "gaussian_k3", c("ID", "BIC")]
  names(baseline)[[2L]] <- "baseline"
  differences <- merge(model_table, baseline, by = "ID")
  differences <- differences[differences$model_key != "gaussian_k3", ]
  differences$delta <- differences$BIC - differences$baseline
  differences$model_key <- factor(
    differences$model_key, levels = c("gaussian_k2", "mixture_k2")
  )
  panel_a <- ggplot2::ggplot(
    differences,
    ggplot2::aes(ID, delta, color = model_key, shape = model_key)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35) +
    ggplot2::geom_point(size = 2.3, position = ggplot2::position_dodge(width = 0.25)) +
    ggplot2::scale_color_manual(values = MODEL_COLORS[c("gaussian_k2", "mixture_k2")], labels = MODEL_LABELS[c("gaussian_k2", "mixture_k2")]) +
    ggplot2::scale_shape_manual(values = c(15, 16), labels = MODEL_LABELS[c("gaussian_k2", "mixture_k2")]) +
    ggplot2::labs(x = "Elk identifier", y = "BIC minus BIC(Gaussian K=3)", title = "Penalized fit") +
    study_theme()

  model_table$model_key <- factor(model_table$model_key, levels = names(MODEL_LABELS))
  panel_b <- ggplot2::ggplot(
    model_table,
    ggplot2::aes(ID, PIT_KS_distance, fill = model_key)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78), width = 0.70) +
    ggplot2::scale_fill_manual(values = MODEL_COLORS, labels = MODEL_LABELS) +
    ggplot2::labs(x = "Elk identifier", y = "One-step PIT distance", title = "Predictive calibration") +
    study_theme()
  combined <- panel_a + panel_b + patchwork::plot_annotation(tag_levels = "A")
  save_elk_plot(combined, "elk_model_diagnostics", 7.4, 3.2)
}

plot_observation_diagnostics <- function(step_data, fits_by_animal) {
  density_rows <- list()
  pit_rows <- list()
  index <- 1L
  for (animal in names(fits_by_animal)) {
    y <- step_data$log1p_step_km[step_data$ID == animal]
    grid_padding <- 0.08 * max(diff(range(y)), 1)
    grid <- seq(min(y) - grid_padding, max(y) + grid_padding, length.out = 500)
    for (model_key in names(MODEL_LABELS)) {
      fit <- fits_by_animal[[animal]][[model_key]]
      density_rows[[index]] <- data.frame(
        ID = animal,
        model_key = model_key,
        value = grid,
        density = occupation_weighted_density(grid, fit)
      )
      pit <- if (model_key == "mixture_k2") pit_mixture(y, fit) else pit_gaussian(y, fit)
      pit_rows[[index]] <- data.frame(ID = animal, model_key = model_key, pit = pit)
      index <- index + 1L
    }
  }
  density_data <- do.call(rbind, density_rows)
  pit_data <- do.call(rbind, pit_rows)
  density_data$model_key <- factor(density_data$model_key, levels = names(MODEL_LABELS))
  pit_data$model_key <- factor(pit_data$model_key, levels = names(MODEL_LABELS))

  panel_a <- ggplot2::ggplot(step_data, ggplot2::aes(log1p_step_km)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 28,
      fill = "grey82",
      color = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_line(
      data = density_data,
      ggplot2::aes(value, density, color = model_key, linetype = model_key),
      linewidth = 0.65,
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(~ID, ncol = 1, scales = "free") +
    ggplot2::scale_color_manual(values = MODEL_COLORS, labels = MODEL_LABELS) +
    ggplot2::scale_linetype_manual(values = c("dashed", "solid", "dotted"), labels = MODEL_LABELS) +
    ggplot2::labs(x = "log(1 + step length in km)", y = "Density", title = "Observed distribution and fitted densities") +
    study_theme(8)

  reference <- data.frame(x = c(0, 1), y = c(0, 1))
  panel_b <- ggplot2::ggplot(
    pit_data[pit_data$model_key != "gaussian_k2", ],
    ggplot2::aes(pit, color = model_key, linetype = model_key)
  ) +
    ggplot2::stat_ecdf(linewidth = 0.65, geom = "step") +
    ggplot2::geom_line(
      data = reference,
      ggplot2::aes(x, y),
      inherit.aes = FALSE,
      color = "grey35",
      linetype = "dashed",
      linewidth = 0.45
    ) +
    ggplot2::facet_wrap(~ID, ncol = 1) +
    ggplot2::scale_color_manual(values = MODEL_COLORS, labels = MODEL_LABELS) +
    ggplot2::scale_linetype_manual(values = c(gaussian_k3 = "solid", mixture_k2 = "dotted"), labels = MODEL_LABELS) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "PIT value", y = "Empirical CDF", title = "Empirical one-step PIT distributions") +
    study_theme(8)
  combined <- panel_a + panel_b + patchwork::plot_annotation(tag_levels = "A")
  save_elk_plot(combined, "elk_observation_diagnostics", 7.4, 8.2)
}

plot_decoding <- function(decoding_table) {
  long <- tidyr::pivot_longer(
    decoding_table,
    cols = c(gaussian_k3_state, mixture_k2_state),
    names_to = "model_key",
    values_to = "state"
  )
  long$model_key <- factor(
    long$model_key,
    levels = c("gaussian_k3_state", "mixture_k2_state"),
    labels = c("Gaussian K=3", "Mixture K=2")
  )
  long$state <- factor(long$state)
  plot <- ggplot2::ggplot(long, ggplot2::aes(step_index, step_km)) +
    ggplot2::geom_line(color = "grey78", linewidth = 0.25) +
    ggplot2::geom_point(ggplot2::aes(color = state), size = 0.65) +
    ggplot2::facet_grid(ID ~ model_key, scales = "free_x") +
    ggplot2::scale_color_manual(values = STATE_COLORS) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.01)) +
    ggplot2::labs(x = "Daily step index", y = "Step length (km)", color = "State") +
    study_theme(8)
  save_elk_plot(plot, "elk_decoding_comparison", 7.4, 7.8)
}

plot_state_mapping <- function(mapping_table) {
  plot <- ggplot2::ggplot(
    mapping_table,
    ggplot2::aes(
      factor(mixture_k2_state, labels = c("aggregate 1", "aggregate 2")),
      factor(gaussian_k3_state, labels = c("G1", "G2", "G3")),
      fill = row_proportion
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", row_proportion)),
      size = 2.6
    ) +
    ggplot2::facet_wrap(~ID, nrow = 1) +
    ggplot2::scale_fill_gradient(low = "white", high = "#2F6B9A", limits = c(0, 1)) +
    ggplot2::labs(
      x = "Enriched-emission state",
      y = "Gaussian K=3 state",
      fill = "Row proportion",
      title = "Mapping of three Gaussian states to two enriched-emission states"
    ) +
    study_theme(8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  save_elk_plot(plot, "elk_state_mapping", 7.4, 2.6)
}

plot_tracks <- function(decoding_table) {
  long <- tidyr::pivot_longer(
    decoding_table,
    cols = c(gaussian_k3_state, mixture_k2_state),
    names_to = "model_key",
    values_to = "state"
  )
  long$model_key <- factor(
    long$model_key,
    levels = c("gaussian_k3_state", "mixture_k2_state"),
    labels = c("Gaussian K=3", "Mixture K=2")
  )
  long$state <- factor(long$state)
  panels <- lapply(unique(interaction(long$ID, long$model_key, drop = TRUE)), function(group) {
    data <- long[interaction(long$ID, long$model_key, drop = TRUE) == group, ]
    x_limits <- range(c(data$x_from, data$x_to))
    y_limits <- range(c(data$y_from, data$y_to))
    span <- max(diff(x_limits), diff(y_limits))
    x_center <- mean(x_limits)
    y_center <- mean(y_limits)
    x_limits <- x_center + c(-0.5, 0.5) * span
    y_limits <- y_center + c(-0.5, 0.5) * span
    ggplot2::ggplot(data) +
      ggplot2::geom_segment(
        ggplot2::aes(x_from, y_from, xend = x_to, yend = y_to, color = state),
        linewidth = 0.55,
        lineend = "round"
      ) +
      ggplot2::scale_color_manual(values = STATE_COLORS, drop = FALSE) +
      ggplot2::coord_equal(xlim = x_limits, ylim = y_limits, expand = FALSE) +
      ggplot2::labs(
        x = NULL,
        y = NULL,
        color = "State",
        title = paste(unique(data$ID), unique(data$model_key), sep = "\n")
      ) +
      study_theme(8) +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 8)
      )
  })
  plot <- patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(title = "Decoded movement paths")
  save_elk_plot(plot, "elk_tracks_decoded", 7.4, 8.2)
}

main <- function() {
  invisible(lapply(
    c(
      OUTPUT_DIRECTORY, FIGURE_DIRECTORY, ARTICLE_FIGURE_DIRECTORY,
      ARTICLE_TABLE_DIRECTORY, SUPPLEMENT_TABLE_DIRECTORY
    ),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
  raw_data <- load_elk_data()
  step_data <- prepare_steps(raw_data)

  fits_by_animal <- list()
  model_rows <- list()
  gaussian_rows <- list()
  mixture_rows <- list()
  mapping_rows_list <- list()
  transition_rows_list <- list()
  decoding_rows <- list()
  animals <- unique(step_data$ID)

  for (animal_index in seq_along(animals)) {
    animal <- animals[[animal_index]]
    data <- step_data[step_data$ID == animal, ]
    y <- data$log1p_step_km
    message(sprintf("Fitting elk %s (%d steps).", animal, length(y)))
    fits <- fit_elk_models(y, animal_index)
    fits_by_animal[[animal]] <- fits

    model_rows[[animal_index]] <- do.call(rbind, lapply(names(fits), function(model_key) {
      model_row(animal, y, model_key, fits[[model_key]])
    }))
    parameters <- parameter_rows(animal, fits)
    gaussian_rows[[animal_index]] <- parameters$gaussian
    mixture_rows[[animal_index]] <- parameters$mixture
    mapping_rows_list[[animal_index]] <- mapping_rows(animal, fits$gaussian_k3, fits$mixture_k2)
    transition_rows_list[[animal_index]] <- transition_rows(animal, fits)

    decoded <- data
    decoded$gaussian_k2_state <- fits$gaussian_k2$viterbi
    decoded$gaussian_k3_state <- fits$gaussian_k3$viterbi
    decoded$mixture_k2_state <- fits$mixture_k2$viterbi
    decoded$mixture_k2_exploratory_probability <- fits$mixture_k2$posterior[, 2L]
    decoding_rows[[animal_index]] <- decoded
  }

  model_table <- do.call(rbind, model_rows)
  model_table <- dplyr::group_by(model_table, ID)
  model_table <- dplyr::mutate(
    model_table,
    delta_BIC_vs_gaussian_k3 = BIC - BIC[model_key == "gaussian_k3"]
  )
  model_table <- dplyr::ungroup(model_table)
  gaussian_table <- do.call(rbind, gaussian_rows)
  mixture_table <- do.call(rbind, mixture_rows)
  mapping_table <- do.call(rbind, mapping_rows_list)
  transition_table <- do.call(rbind, transition_rows_list)
  decoding_table <- do.call(rbind, decoding_rows)

  utils::write.csv(model_table, file.path(OUTPUT_DIRECTORY, "elk_model_comparison.csv"), row.names = FALSE)
  utils::write.csv(gaussian_table, file.path(OUTPUT_DIRECTORY, "elk_gaussian_state_parameters.csv"), row.names = FALSE)
  utils::write.csv(mixture_table, file.path(OUTPUT_DIRECTORY, "elk_mixture_state_parameters.csv"), row.names = FALSE)
  utils::write.csv(mapping_table, file.path(OUTPUT_DIRECTORY, "elk_state_mapping.csv"), row.names = FALSE)
  utils::write.csv(transition_table, file.path(OUTPUT_DIRECTORY, "elk_transition_matrices.csv"), row.names = FALSE)
  archived_decoding_table <- dplyr::select(
    decoding_table,
    -x_from, -y_from, -x_to, -y_to
  )
  utils::write.csv(
    archived_decoding_table,
    file.path(OUTPUT_DIRECTORY, "elk_decoding.csv"),
    row.names = FALSE
  )

  metadata <- list(
    implementation = "R",
    r_version = R.version.string,
    seed = SEED,
    elk_data_source = "moveHMM::elk_data",
    moveHMM_version = as.character(utils::packageVersion("moveHMM")),
    elk_data_object_sha256 = EXPECTED_DATA_DIGEST,
    raw_data_csv_distributed = FALSE,
    response = "log(1 + daily step length in kilometres)",
    n_starts = N_STARTS,
    screen_iterations = SCREEN_ITER,
    refined_starts = REFINE_TOP,
    maximum_refinement_iterations = MAX_ITER,
    minimum_standard_deviation = MINIMUM_SD,
    animals = animals,
    n_steps = stats::setNames(as.list(as.integer(table(step_data$ID)[animals])), animals),
    models = as.list(MODEL_LABELS),
    enriched_component_counts = c(2L, 2L)
  )
  jsonlite::write_json(
    metadata,
    file.path(OUTPUT_DIRECTORY, "elk_analysis_metadata.json"),
    auto_unbox = TRUE,
    digits = 16,
    pretty = TRUE
  )

  write_model_table(model_table)
  write_parameter_table(mixture_table)
  write_predictive_table(model_table)
  plot_model_diagnostics(model_table)
  plot_observation_diagnostics(step_data, fits_by_animal)
  plot_decoding(decoding_table)
  plot_state_mapping(mapping_table)
  plot_tracks(decoding_table)

  print(model_table[c(
    "ID", "model", "logLik", "n_parameters", "BIC",
    "delta_BIC_vs_gaussian_k3", "PIT_KS_distance", "PIT_90pct_coverage"
  )])
}

if (sys.nframe() == 0L) {
  main()
}
