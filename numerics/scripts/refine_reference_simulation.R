library(Rcpp);library(ggplot2);library(dplyr);library(jsonlite)
source('numerics/R/hmm_em.R');source('numerics/R/mixture_hmm_em.R');source('numerics/R/simulate.R')
source('numerics/R/diagnostics.R');source('numerics/R/plotting.R');source('numerics/R/optimization_audit.R')
simulation<-simulate_hmm(1500L,20260610L);y<-simulation$y
fs<-fit_simulation_audited(y,20260610L)
# Public downstream plot/diagnostic interfaces expect vectors for Gaussian emissions.
as_gaussian<-function(f){
 o<-order(f$means[,1]);f$initial<-f$initial[o];f$transition<-f$transition[o,o];f$means<-f$means[o,1];f$sds<-f$sds[o,1];f$posterior<-f$posterior[,o]
 f$viterbi<-viterbi_gaussian(y,f$initial,f$transition,f$means,f$sds);class(f)<-c('gaussian_hmm_fit','hmm_fit');f
}
g2<-as_gaussian(fs$G2);g3<-as_gaussian(fs$G3);m<-fs$M22
ord<-order(mixture_state_means(m));m$initial<-m$initial[ord];m$transition<-m$transition[ord,ord];m$weights<-m$weights[ord,];m$means<-m$means[ord,];m$sds<-m$sds[ord,];m$component_counts<-m$component_counts[ord];m$posterior<-m$posterior[,ord];m$viterbi<-viterbi_mixture(y,m$initial,m$transition,m$weights,m$means,m$sds);class(m)<-c('mixture_hmm_fit','hmm_fit')
saveRDS(list(G2=g2,G3=g3,M22=m), 'numerics/output/optimization_audit/reference_fits.rds')
p<-parse('numerics/scripts/run_review_simulation.R');eval(p[-length(p)])
comparison<-model_rows(g2,g3,m,1500);write_csv(comparison,'numerics/output/review_model_comparison.csv');write_article_model_table(comparison)
write_csv(rbind(gaussian_fit_summary(g2,1500),gaussian_fit_summary(g3,1500)),'numerics/output/review_gaussian_fit_summary.csv')
write_csv(mixture_fit_summary(m,1500),'numerics/output/review_mixture_fit_summary.csv')
plot_simulation_density(y,g2,g3,m,'article/figures')
t2<-crosstab_true_fitted(simulation$refined_state,g2);t3<-crosstab_true_fitted(simulation$refined_state,g3)
write_csv(as.data.frame.matrix(t2),'numerics/output/crosstab_fit2.csv',row.names=TRUE);write_csv(as.data.frame.matrix(t3),'numerics/output/crosstab_fit3.csv',row.names=TRUE)
plot_decoding_heatmap(t2,t3,'article/figures')
plot_corrective_posterior(y,g3,which.min(abs(g3$means-TRUE_REFINED_MEANS[2])),'article/figures')
a2<-residual_acf(one_step_residuals_gaussian(y,g2));a3<-residual_acf(one_step_residuals_gaussian(y,g3))
write_csv(a2,'numerics/output/residual_acf_fit2.csv');write_csv(a3,'numerics/output/residual_acf_fit3.csv');plot_residual_acf(a2,a3,'article/figures')
print(comparison)
