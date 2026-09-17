# ggplot2 figures for the numerical study.

save_ggplot <- function(plot, directory, name, width, height) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file.path(directory, paste0(name, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    filename = file.path(directory, paste0(name, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
  invisible(plot)
}

study_theme <- function(base_size = 9) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      plot.title.position = "plot"
    )
}

occupation_weighted_density <- function(grid, fit) {
  occupation <- colMeans(fit$posterior)
  if (inherits(fit, "mixture_hmm_fit")) {
    output <- numeric(length(grid))
    for (state in seq_len(fit$K)) {
      for (component in seq_len(fit$component_counts[[state]])) {
        output <- output + occupation[[state]] * fit$weights[state, component] *
          stats::dnorm(
            grid,
            mean = fit$means[state, component],
            sd = fit$sds[state, component]
          )
      }
    }
    return(output)
  }
  Reduce(
    `+`,
    lapply(seq_len(fit$K), function(state) {
      occupation[[state]] * stats::dnorm(
        grid, mean = fit$means[[state]], sd = fit$sds[[state]]
      )
    })
  )
}

plot_simulation_density <- function(y, fit2, fit3, mixture2, figure_directory) {
  grid <- seq(min(y) - 0.8, max(y) + 0.8, length.out = 600)
  density_data <- rbind(
    data.frame(
      observation = grid,
      density = occupation_weighted_density(grid, fit2),
      model = "Gaussian HMM, K=2"
    ),
    data.frame(
      observation = grid,
      density = occupation_weighted_density(grid, fit3),
      model = "Gaussian HMM, K=3"
    ),
    data.frame(
      observation = grid,
      density = occupation_weighted_density(grid, mixture2),
      model = "Mixture-emission HMM, K=2"
    )
  )
  density_data$model <- factor(
    density_data$model,
    levels = c(
      "Gaussian HMM, K=2",
      "Gaussian HMM, K=3",
      "Mixture-emission HMM, K=2"
    )
  )
  plot <- ggplot2::ggplot(data.frame(observation = y), ggplot2::aes(observation)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 60,
      fill = "grey70",
      color = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_line(
      data = density_data,
      ggplot2::aes(y = density, color = model, linetype = model),
      linewidth = 0.8
    ) +
    ggplot2::scale_color_manual(values = c("#C76E2E", "#2F6B9A", "#2D8A66")) +
    ggplot2::scale_linetype_manual(values = c("dashed", "solid", "dotted")) +
    ggplot2::labs(x = "Observation", y = "Density") +
    study_theme()
  save_ggplot(plot, figure_directory, "simulation_density", 6.2, 3.7)
}

plot_decoding_heatmap <- function(table2, table3, figure_directory) {
  make_data <- function(tab, K) {
    proportions <- prop.table(tab, margin = 1)
    output <- as.data.frame(proportions, responseName = "proportion")
    names(output)[1:2] <- c("true_refined", "fitted_state")
    output$K <- paste0("Gaussian HMM, K=", K)
    output
  }
  data <- rbind(make_data(table2, 2), make_data(table3, 3))
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(fitted_state, true_refined, fill = proportion)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", proportion)),
      size = 2.8
    ) +
    ggplot2::facet_wrap(~K, scales = "free_x") +
    ggplot2::scale_fill_gradient(low = "white", high = "#2F6B9A", limits = c(0, 1)) +
    ggplot2::labs(
      x = "Fitted state",
      y = "True refined state",
      fill = "Row proportion"
    ) +
    study_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  save_ggplot(plot, figure_directory, "state_decoding_heatmap", 7.2, 3.0)
}

plot_corrective_posterior <- function(
    y, fit3, corrective_state, figure_directory) {
  breaks <- unique(as.numeric(stats::quantile(
    y, seq(0, 1, length.out = 41), names = FALSE
  )))
  bin <- cut(y, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  point_data <- data.frame(
    observation = y,
    posterior = fit3$posterior[, corrective_state],
    bin = bin
  )
  binned <- dplyr::summarise(
    dplyr::group_by(point_data, bin),
    observation = mean(observation),
    posterior = mean(posterior),
    n = dplyr::n(),
    .groups = "drop"
  )
  binned <- dplyr::filter(binned, n >= 5)
  plot <- ggplot2::ggplot(
    point_data,
    ggplot2::aes(observation, posterior)
  ) +
    ggplot2::geom_point(color = "#2F6B9A", alpha = 0.30, size = 0.8) +
    ggplot2::geom_line(
      data = binned,
      color = "#C44E52",
      linewidth = 0.8
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = expression("Observation " * Y[t]),
      y = "Posterior probability of corrective state"
    ) +
    study_theme()
  save_ggplot(plot, figure_directory, "posterior_corrective_state", 6.2, 3.5)
}

plot_bic_boxplot <- function(delta_bic, figure_directory) {
  frequency <- mean(delta_bic < 0)
  plot <- ggplot2::ggplot(
    data.frame(group = "replications", delta_bic = delta_bic),
    ggplot2::aes(group, delta_bic)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "#C44E52", linetype = "dashed") +
    ggplot2::geom_boxplot(width = 0.4, fill = "#D8E6F2", color = "#2F6B9A") +
    ggplot2::labs(
      x = NULL,
      y = expression(BIC(K == 3) - BIC(K == 2)),
      title = sprintf("BIC favours K=3 in %.1f%% of replications", 100 * frequency)
    ) +
    study_theme()
  save_ggplot(plot, figure_directory, "bic_boxplot", 4.4, 3.5)
}

plot_bic_by_sample_size <- function(summary, figure_directory) {
  data <- rbind(
    data.frame(
      T = summary$T,
      selection_frequency = summary$prop_gaussian_BIC_favours_K3,
      comparison = "Gaussian: BIC favours K=3"
    ),
    data.frame(
      T = summary$T,
      selection_frequency = summary$prop_enriched_BIC_favours_K2_mix,
      comparison = "Enriched: BIC favours K=2 mixture"
    )
  )
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(factor(T), selection_frequency, fill = comparison)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.68) +
    ggplot2::scale_fill_manual(values = c("#2D8A66", "#2F6B9A")) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "T", y = "Selection frequency") +
    study_theme()
  save_ggplot(plot, figure_directory, "bic_by_T", 5.2, 3.4)
}

plot_residual_acf <- function(acf2, acf3, figure_directory) {
  data <- rbind(
    transform(acf2, model = "Gaussian HMM, K=2"),
    transform(acf3, model = "Gaussian HMM, K=3")
  )
  plot <- ggplot2::ggplot(data, ggplot2::aes(lag, acf, color = model)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_color_manual(values = c("#C76E2E", "#2F6B9A")) +
    ggplot2::labs(x = "Lag", y = "Residual ACF") +
    study_theme()
  save_ggplot(plot, figure_directory, "residual_acf", 5.4, 3.4)
}
