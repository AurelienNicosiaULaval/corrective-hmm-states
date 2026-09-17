# A domain-informed comparator for unphased, heterozygous SNP measurements.
# Each state has a symmetric pair of normal densities around BAF = 1/2.
# On |BAF - 1/2| this is a folded-normal HMM. Scores returned here are on BAF.
library(Rcpp)

sym_simplex <- function(z,lower=1e-6) {
  free<-rep(TRUE,length(z));p<-rep(lower,length(z))
  if(sum(z)<=0)return(rep(1/length(z),length(z)))
  repeat {
    p[free]<-z[free]/sum(z[free])*(1-sum(!free)*lower)
    bad<-free & p<lower
    if(!any(bad))break
    free[bad]<-FALSE;p[bad]<-lower
  }
  p
}
sym_evaluate <- function(y,start,f,decode=FALSE) {
  x<-y-.5;k<-length(f$mu)
  plus<-vapply(seq_len(k),function(j)dnorm(x,f$mu[j],f$sd[j]),numeric(length(x)))
  minus<-vapply(seq_len(k),function(j)dnorm(x,-f$mu[j],f$sd[j]),numeric(length(x)))
  emission<-pmax((plus+minus)/2,1e-300)
  out<-fb_sequences_cpp(emission,f$initial,f$transition,start)
  out$positive_responsibility<-out$posterior*plus/(2*emission)
  out$negative_responsibility<-out$posterior*minus/(2*emission)
  out$log_scores<-log(out$scales)
  if(decode)out$viterbi<-viterbi_sequences_cpp(log(emission),f$initial,f$transition,start)
  out
}
sym_fit <- function(y,start,f,max_iter=1500,tolerance=1e-8,floor=.01) {
  trace<-numeric();converged<-FALSE;x<-y-.5
  for(iter in seq_len(max_iter)) {
    e<-sym_evaluate(y,start,f);trace<-c(trace,e$log_likelihood)
    if(length(trace)>1) {
      delta<-tail(diff(trace),1)
      if(delta< -1e-7*(1+abs(e$log_likelihood)))stop('Symmetric EM likelihood decreased')
      if(delta>=-1e-8 && delta<tolerance*(1+abs(e$log_likelihood))){converged<-TRUE;break}
    }
    f$initial<-sym_simplex(e$initial_sum)
    f$transition<-t(apply(e$transition_sum,1,sym_simplex))
    a<-e$positive_responsibility;b<-e$negative_responsibility
    mass<-colSums(a+b)
    signed_mean<-colSums((a-b)*x)/mass
    for(j in seq_along(f$mu)) {
      f$sd[j]<-sqrt(max(floor^2,sum(a[,j]*(x-signed_mean[j])^2+b[,j]*(x+signed_mean[j])^2)/mass[j]))
    }
    f$mu<-abs(signed_mean)
  }
  e<-sym_evaluate(y,start,f,TRUE)
  f$log_likelihood<-e$log_likelihood;f$posterior<-e$posterior;f$viterbi<-e$viterbi
  f$trace<-c(trace,e$log_likelihood);f$converged<-converged
  k<-length(f$mu);f$n_parameters<-k*k-1+2*k
  f$BIC<- -2*f$log_likelihood+f$n_parameters*log(length(y));f
}
