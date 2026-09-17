library(Rcpp)
source('numerics/R/hmm_em.R'); source('numerics/R/mixture_hmm_em.R')
source('numerics/R/fast_em.R')
set.seed(32); y <- rnorm(73)
for(counts in list(c(1L,1L),c(2L,1L),c(2L,2L),c(2L,2L,1L))) {
 for(lengths in list(NULL,c(31L,42L))) {
  f <- initial_mixture_parameters(y,length(counts),max(counts),counts,.05)
  for(it in c(1L,5L,25L)) {
   r <- mixture_em_loop(y,f$initial,f$transition,f$weights,f$means,f$sds,counts,it,1e-9,.05,lengths)
   c <- fast_em(y,f,counts,lengths,it)
   stopifnot(max(abs(r$log_likelihood_trace-c$log_likelihood_trace))<1e-8,
             max(abs(r$posterior-c$posterior))<1e-8,
             max(abs(r$transition-c$transition))<1e-8,
             max(abs(r$means-c$means))<1e-8)
  }
 }
}
cat('Compiled EM matches R for four models, both boundary settings and three iteration counts.\n')
