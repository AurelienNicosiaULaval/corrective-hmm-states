library(Rcpp)
sourceCpp('genomics/scripts/gaussian_mixture_hmm.cpp')
source('genomics/scripts/symmetric_genomic_hmm.R')
y<-c(.19,.72,.48,.53,.82,.26,.43,.61)
start<-c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE)
f<-list(mu=c(.02,.28),sd=c(.06,.09),initial=c(.6,.4),
        transition=matrix(c(.91,.09,.13,.87),2,byrow=TRUE))
a<-sym_evaluate(y,start,f,TRUE)
# Allele designation can be exchanged independently at any subset of markers.
flip<-c(TRUE,FALSE,TRUE,TRUE,FALSE,TRUE,FALSE,TRUE)
y2<-y;y2[flip]<-1-y2[flip]
b<-sym_evaluate(y2,start,f,TRUE)
stopifnot(abs(a$log_likelihood-b$log_likelihood)<1e-12,
          max(abs(a$posterior-b$posterior))<1e-12,identical(a$viterbi,b$viterbi))
# Direct path enumeration is independent of the forward-backward code.
direct<-function(z) {
  paths<-as.matrix(expand.grid(rep(list(1:2),length(z))))
  probs<-apply(paths,1,function(s) {
    density<-(dnorm(z-.5,f$mu[s],f$sd[s])+dnorm(z-.5,-f$mu[s],f$sd[s]))/2
    f$initial[s[1]]*prod(f$transition[cbind(head(s,-1),tail(s,-1))])*prod(density)
  })
  log(sum(probs))
}
stopifnot(abs(a$log_likelihood-direct(y[1:4])-direct(y[5:8]))<1e-12)
fit<-sym_fit(y,start,f)
stopifnot(all(diff(fit$trace)>-1e-8),all(fit$sd>=.01))
cat('Passed: direct path likelihood, multiple sequence resets, arbitrary allele-flip invariance, posterior and Viterbi invariance, and EM ascent.\n')
