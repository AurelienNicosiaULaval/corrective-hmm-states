#!/usr/bin/env Rscript
# Rebuild manuscript and supplementary tables from saved numerical results.
library(data.table)

write_table <- function(path, caption, label, columns, heading, rows) {
  writeLines(c('\\begin{table}[htbp]', '\\centering', '\\small',
               paste0('\\caption{', caption, '}'), paste0('\\label{', label, '}'),
               paste0('\\begin{tabular}{', columns, '}'), '\\toprule',
               paste0(heading, ' \\\\'), '\\midrule', rows,
               '\\bottomrule', '\\end{tabular}', '\\end{table}'), path)
}
all_models <- c('G2', 'G3', 'G4', 'G5', 'M21', 'M221', 'F2', 'F3')
metrics <- rbindlist(list(fread('genomics/results/genomic_pilot/metrics.csv'),
                        fread('genomics/results/genomic_symmetry_benchmark/metrics.csv')), fill = TRUE)
stopifnot(!anyDuplicated(metrics$model), setequal(metrics$model, all_models))
metrics <- metrics[match(all_models, model)]
scores <- fread('genomics/results/genomic_validation/scores.csv')
score <- function(model_name, dataset_name) {
  result <- scores[model == model_name & dataset == dataset_name, log_score]
  stopifnot(length(result) == 1L, is.finite(result))
  result
}
main <- metrics[model %in% c('G2', 'G3', 'G4', 'G5', 'M21', 'M221')]
rows <- vapply(seq_len(nrow(main)), function(i) {
  z <- main[i]
  sprintf('%s & %d & %d & %d & %.1f & %.2f \\\\', z$model, z$states,
          z$components, z$parameters, z$BIC, score(z$model, 'replication'))
}, character(1))
write_table('article/tables/genomic_comparison.tex',
  'Chromosome-1 comparison. $K$ is the number of Markov states, $R$ the total number of Gaussian components, and $p$ the number of free parameters. Lower BIC is preferred. Repeat scores are sequential log scores on chromosome 1 of the second array pair, with fitted parameters held fixed; higher scores are preferred. This is a technical repeat of the same cell lines.',
  'tab:genomic-comparison', 'lrrrrr',
  'Model & $K$ & $R$ & $p$ & BIC & Repeat score', rows)

rows <- vapply(seq_len(nrow(metrics)), function(i) {
  z <- metrics[i]
  sprintf('%s & %d & %.1f & %.2f & %.2f & %.2f \\\\', z$model,
          z$parameters, z$BIC, score(z$model, 'replication'),
          score(z$model, 'validation'), score(z$model, 'confirmation'))
}, character(1))
write_table('supplement/tables/genomic_comparison_full.tex',
  'Full genomic comparison. All parameters are estimated on chromosome 1. Lower BIC and higher sequential log scores are preferred in their respective columns. F2/F3 impose reflection symmetry. G4 was added retrospectively after the original seven-model evaluation; it was not part of the original model freeze.',
  'tab:si-genomic-comparison', 'lrrrrr',
  'Model & $p$ & BIC & Repeat & Chr. 12--17 & Chr. 18--22', rows)

rows <- vapply(seq_len(nrow(metrics)), function(i) {
  z <- metrics[i]
  sprintf('%s & %.3f & %d & %s & %d \\\\', z$model, z$log_likelihood,
          z$runs, ifelse(z$converged, 'Yes', 'No'), z$floor_hits)
}, character(1))
write_table('supplement/tables/genomic_estimation.tex',
  'Selected estimates on chromosome 1. Runs count consecutive Viterbi labels within each sequence; they are not verified genomic alterations. Convergence refers to the final relative likelihood stopping criterion, not to proof of a global optimum.',
  'tab:si-genomic-estimation', 'lrrrr',
  'Model & $\\ell$ & Runs & Converged & Floor hits', rows)

by_chr <- fread('genomics/results/genomic_validation/chromosome_scores.csv')
wide <- dcast(by_chr[dataset != 'replication'], dataset + chromosome + n ~ model,
              value.var = 'log_score')
setorder(wide, chromosome)
rows <- vapply(seq_len(nrow(wide)), function(i) {
  z <- wide[i]
  sprintf('%d & %d & %.2f & %.2f & %.2f \\\\', z$chromosome, z$n,
          z$G4 - z$M221, z$G5 - z$M221, z$F3 - z$M221)
}, character(1))
write_table('supplement/tables/genomic_chromosomes.tex',
  'Sequential log-score differences by excluded chromosome, using fixed chromosome-1 estimates. Differences are comparator minus M221, so positive values favour the comparator. The G4 comparison was added retrospectively.',
  'tab:si-genomic-chromosomes', 'rrrrr',
  'Chromosome & $n$ & G4 $-$ M221 & G5 $-$ M221 & F3 $-$ M221', rows)

calibration <- fread('genomics/results/genomic_diagnostics/predictive_calibration.csv')
calibration <- calibration[order(match(dataset, c('validation', 'confirmation')),
                                 match(model, all_models))]
rows <- vapply(seq_len(nrow(calibration)), function(i) {
  z <- calibration[i]
  sprintf('%s & %s & %.2f & %.3f & %.3f & %.3f \\\\',
          ifelse(z$dataset == 'validation', '12--17', '18--22'), z$model,
          100 * z$coverage95, z$pit_KS_distance, z$residual_lag1,
          z$absolute_residual_lag1)
}, character(1))
write_table('supplement/tables/genomic_calibration.tex',
  'Fixed-parameter predictive diagnostics. Coverage is the percentage inside the central 95\\% predictive interval; $D_{\\mathrm{PIT}}$ is the empirical distance from a uniform PIT distribution. The last two columns are lag-one correlations of normal-score residuals and their absolute values, excluding sequence boundaries.',
  'tab:si-genomic-calibration', 'llrrrr',
  'Chromosomes & Model & Coverage (\\%) & $D_{\\mathrm{PIT}}$ & $r_1$ & $r_1^{|\\cdot|}$', rows)

sensitivity <- fread('genomics/results/genomic_sensitivity/summary.csv')
labels <- c(floor_005 = '$\\sigma_{\\min}=0.005$',
            floor_020 = '$\\sigma_{\\min}=0.020$', thin_2 = 'Every second marker',
            thin_5 = 'Every fifth marker', normal_025 = 'Normal fraction $[0.25,0.75]$',
            normal_035 = 'Normal fraction $[0.35,0.65]$')
rows <- vapply(seq_len(nrow(sensitivity)), function(i) {
  z <- sensitivity[i]
  sprintf('%s & %d & %.2f & %d \\\\', labels[z$setting], z$n,
          100 * z$state_agreement, z$runs)
}, character(1))
write_table('supplement/tables/genomic_sensitivity.tex',
  'M221 sensitivity. Agreement is evaluated at markers shared with the primary analysis; parameters are refitted within each scenario.',
  'tab:si-genomic-sensitivity', 'lrrr',
  'Scenario & $n$ & Agreement (\\%) & Runs', rows)
cat('Six publication tables regenerated from saved results.\n')
