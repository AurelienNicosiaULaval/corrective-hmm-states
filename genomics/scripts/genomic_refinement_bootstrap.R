# Parametric bootstrap for the prespecified (2,2,1) refinement constraints.
# Conditional on the fitted null model and the two observed sequence lengths.
library(Rcpp); library(data.table)
source('numerics/R/hmm_em.R');source('numerics/R/mixture_hmm_em.R')
source('numerics/R/fast_em.R');source('numerics/R/revision_validation.R')
out<-'genomics/results/refinement_bootstrap';dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-fread('genomics/data/genomics/development_chr1.csv');lens<-as.integer(rle(d$sequence)$lengths)
null0<-readRDS('genomics/results/genomic_pilot/M221.rds');alt0<-readRDS('genomics/results/genomic_pilot/G5.rds')
finish<-function(y,f,counts){
 for(i in 1:8){f<-fast_em(y,f,counts,lens,3000L,1e-8,.01);if(relative_gain(f)<=1e-8)break}
 if(relative_gain(f)>1.01e-8)stop('Iteration cap before convergence')
 f
}
fit_pair<-function(y,seed){
 set.seed(seed)
 fit<-function(counts,extra){
  candidates<-c(extra,lapply(1:10,function(i)initial_mixture_parameters(y,length(counts),max(counts),counts,.01)))
  screened<-lapply(candidates,function(f)fast_em(y,f,counts,lens,100L,1e-5,.01))
  ids<-head(order(vapply(screened,`[[`,numeric(1),'log_likelihood'),decreasing=TRUE),4)
  refined<-lapply(screened[ids],function(f)finish(y,f,counts))
  refined[[which.max(vapply(refined,`[[`,numeric(1),'log_likelihood'))]]
 }
 n<-fit(c(2L,2L,1L),list(null0))
 expanded<-expand_mixture_fit(n)
 a<-fit(rep(1L,5L),list(expanded,alt0))
 # A feasible exact refinement is never discarded at the screening stage.
 check<-finish(y,expanded,rep(1L,5L));if(check$log_likelihood>a$log_likelihood)a<-check
 stopifnot(a$log_likelihood>=n$log_likelihood-1e-6)
 list(null=n,alternative=a,statistic=max(0,2*(a$log_likelihood-n$log_likelihood)))
}
obsfile<-file.path(out,'observed.rds')
if(file.exists(obsfile))obs<-readRDS(obsfile) else {
 obs<-fit_pair(d$beta_corrected,202609170L)
 saveRDS(obs,obsfile)
}
print(c(observed_statistic=obs$statistic,null_gain=obs$null$log_likelihood-null0$log_likelihood,alternative_gain=obs$alternative$log_likelihood-alt0$log_likelihood))
generator<-obs$null
simulate_null<-function(seed){
 set.seed(seed);ys<-lapply(lens,function(n){
  st<-integer(n);st[1]<-sample.int(3,1,prob=generator$initial)
  for(t in 2:n)st[t]<-sample.int(3,1,prob=generator$transition[st[t-1],])
  comp<-vapply(st,function(k)sample.int(generator$component_counts[k],1,prob=generator$weights[k,seq_len(generator$component_counts[k])]),integer(1))
  ix<-cbind(st,comp);rnorm(n,generator$means[ix],generator$sds[ix])
 });unlist(ys,use.names=FALSE)
}
B<-as.integer(Sys.getenv('BOOTSTRAP_B','199'))
run<-function(b){
 path<-file.path(out,sprintf('replicate_%03d.rds',b));if(file.exists(path))return(path)
 seed<-202609170L+b;y<-simulate_null(seed);fit<-fit_pair(y,seed+10000L)
 for(n in c('null','alternative'))fit[[n]]$posterior<-NULL
 fit$seed<-seed;saveRDS(fit,paste0(path,'.tmp'));file.rename(paste0(path,'.tmp'),path)
 cat('Bootstrap ',b,'/',B,': ',fit$statistic,'\n',sep='');path
}
paths<-parallel::mclapply(seq_len(B),run,mc.cores=as.integer(Sys.getenv('BOOTSTRAP_WORKERS','4')),mc.preschedule=FALSE)
if(any(vapply(paths,inherits,logical(1),'try-error')))stop('Bootstrap fit failure')
x<-vapply(paths,function(p)readRDS(p)$statistic,numeric(1));exceed<-sum(x>=obs$statistic)
p<-(1+exceed)/(B+1);ci<-binom.test(exceed,B)$conf.int
summary<-data.frame(statistic=obs$statistic,B=B,exceedances=exceed,p_bootstrap=p,monte_carlo_se=sqrt(p*(1-p)/(B+1)),binomial_lower=ci[1],binomial_upper=ci[2],null_log_likelihood=obs$null$log_likelihood,alternative_log_likelihood=obs$alternative$log_likelihood)
write.csv(data.frame(replication=seq_len(B),statistic=x),file.path(out,'replications.csv'),row.names=FALSE)
write.csv(summary,file.path(out,'summary.csv'),row.names=FALSE);print(summary)
writeLines(c('Model-conditional parametric bootstrap; B=199, plus-one tail probability.',
 'Two independent sequences retain their observed lengths and a shared freely estimated initial distribution.',
 'Each fit screens 10 random starts plus saved-data starts for 100 iterations; best four refined at 1e-8.',
 'The G5 fit additionally retains/refines the exact M221 expansion. Fixed state grouping: (2,2,1).',
 'This calibration assumes the fitted Gaussian-mixture HMM, not independent identically distributed markers.',
 'Residual dependence and other misspecification can invalidate model-conditional calibration; nonrejection is not equivalence.'),file.path(out,'protocol.txt'))
