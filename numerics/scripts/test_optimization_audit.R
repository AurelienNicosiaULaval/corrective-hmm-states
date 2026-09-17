library(testthat);library(Rcpp)
source('numerics/R/hmm_em.R');source('numerics/R/mixture_hmm_em.R');source('numerics/R/simulate.R');source('numerics/R/revision_validation.R')
sourceCpp('numerics/R/forward_backward.cpp');forward_backward_core<-forward_backward_compiled
r<-read.csv('numerics/output/review_replications.csv');old<-read.csv('numerics/output/optimization_audit/original_replications.csv')
test_that('all 600 fits improve the old likelihoods and satisfy nesting and stopping rules',{
 expect_equal(nrow(r),600L);expect_equal(table(r$T),table(old$T));expect_false(anyDuplicated(r[c('T','rep')])>0)
 for(i in seq_len(nrow(r))){
  row<-r[i,];z<-readRDS(sprintf('numerics/output/optimization_audit/fit_%d_%03d.rds',row$T,row$rep))
  f<-z$fits;y<-simulate_hmm(row$T,row$seed)$y
  for(n in names(f)){
   fit<-f[[n]];ll<-forward_backward_mixture(y,fit$initial,fit$transition,fit$weights,fit$means,fit$sds)$log_likelihood
   expect_equal(ll,fit$log_likelihood,tolerance=1e-10)
   expect_lte(relative_gain(fit),1.01e-9)
   expect_true(all(diff(fit$log_likelihood_trace)> -1e-6))
  }
  expect_gte(f$G3$log_likelihood,f$M21$log_likelihood-1e-6)
  expect_gte(f$M22$log_likelihood,f$M21$log_likelihood-1e-6)
  baseline<-old[old$T==row$T & old$rep==row$rep,]
  for(pair in list(c('G2','logLik_gauss_K2'),c('G3','logLik_gauss_K3'),c('M22','logLik_mix_K2_M2'))){
   expect_equal(f[[pair[1]]]$log_likelihood,row[[pair[2]]],tolerance=1e-10)
   expect_gte(row[[pair[2]]],baseline[[pair[2]]]-1e-5)
  }
  expect_equal(row$BIC_gauss_K2,fit_bic(f$G2,row$T),tolerance=1e-10)
  expect_equal(row$BIC_gauss_K3,fit_bic(f$G3,row$T),tolerance=1e-10)
  expect_equal(row$BIC_mix_K2_M2,fit_bic(f$M22,row$T),tolerance=1e-10)
 }
})
test_that('summary frequencies agree with row-level comparisons',{
 s<-read.csv('numerics/output/review_summary_by_T.csv')
 for(i in seq_len(nrow(s))){x<-r[r$T==s$T[i],];expect_equal(s$prop_gaussian_BIC_favours_K3[i],mean(x$BIC_gauss_K3<x$BIC_gauss_K2));expect_equal(s$prop_enriched_BIC_favours_K2_mix[i],mean(x$BIC_mix_K2_M2<x$BIC_gauss_K3))}
})
if(file.exists('genomics/results/refinement_bootstrap/summary.csv'))test_that('bootstrap preserves nesting, convergence and its documented tail calculation',{
 p<-'genomics/results/refinement_bootstrap';s<-read.csv(file.path(p,'summary.csv'));rs<-read.csv(file.path(p,'replications.csv'))
 expect_equal(nrow(rs),199L);expect_equal(s$p_bootstrap,(1+sum(rs$statistic>=s$statistic))/200)
 for(i in 1:199){f<-readRDS(file.path(p,sprintf('replicate_%03d.rds',i)));expect_gte(f$alternative$log_likelihood,f$null$log_likelihood-1e-6);expect_lte(relative_gain(f$alternative),1.01e-8);expect_lte(relative_gain(f$null),1.01e-8);expect_equal(f$statistic,rs$statistic[i],tolerance=1e-10)}
})
cat('Corrected simulation and available bootstrap outputs verified.\n')
