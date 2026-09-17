# Data-derived starts for the simulation. No generating parameter is used here.
source('numerics/R/fast_em.R')
source('numerics/R/revision_validation.R')
best_fit <- function(fs) fs[[which.max(vapply(fs,`[[`,numeric(1),'log_likelihood'))]]
# Preserve the original Gaussian random-start distribution and RNG sequence.
gaussian_em_loop <- function(y, initial, transition, means, sds, max_iter,
                             tolerance, sigma_min, lengths=NULL) {
  fast_em(y,list(initial=initial,transition=transition,weights=matrix(1,length(initial),1),
                means=matrix(means,ncol=1),sds=matrix(sds,ncol=1)),
          rep(1L,length(initial)),lengths,max_iter,tolerance,sigma_min)
}
random_fit <- function(y,counts,seed,gaussian=FALSE) {
 set.seed(seed)
 screened <- lapply(seq_len(10L),function(i) {
  if(gaussian) gaussian_em_start(y,length(counts),25L,1e-4,.05) else
   fast_em(y,initial_mixture_parameters(y,length(counts),max(counts),counts,.05),
           counts,max_iter=30L,tolerance=1e-4)
 })
 ids <- head(order(vapply(screened,`[[`,numeric(1),'log_likelihood'),decreasing=TRUE),3L)
 best_fit(lapply(screened[ids],function(f) fast_em(y,f)))
}
split_gaussian <- function(g,state) {
 f<-g; k<-g$K;f$weights<-matrix(0,k,2); f$weights[,1]<-1
 f$means<-cbind(g$means[,1],g$means[,1]); f$sds<-cbind(g$sds[,1],g$sds[,1])
 f$component_counts<-rep(1L,k);f$component_counts[state]<-2L
 f$weights[state,]<-.5
 f$means[state,]<-g$means[state,1]+c(-.7,.7)*g$sds[state,1]
 f$sds[state,]<-g$sds[state,1]*sqrt(1-.7^2);f$M<-2L;f
}
duplicate_singleton <- function(f) {
 s<-which(f$component_counts==1L); f$component_counts[s]<-2L
 f$weights[s,]<-.5;f$means[s,2]<-f$means[s,1];f$sds[s,2]<-f$sds[s,1];f
}
reduce_m22 <- function(f,s) {
 mean<-sum(f$weights[s,]*f$means[s,]);v<-sum(f$weights[s,]*(f$sds[s,]^2+f$means[s,]^2))-mean^2
 f$weights[s,]<-c(1,0);f$means[s,1]<-mean;f$sds[s,1]<-sqrt(max(.05^2,v));f$component_counts[s]<-1L;f
}
fit_simulation_audited <- function(y,seed) {
 g2<-random_fit(y,c(1L,1L),seed+2L,TRUE)
 g3<-random_fit(y,rep(1L,3),seed+3L,TRUE)
 m22<-random_fit(y,c(2L,2L),seed+20L)
 m21<-best_fit(c(list(random_fit(y,c(2L,1L),seed+21L)),
                 lapply(1:2,function(s)fast_em(y,split_gaussian(g2,s))),
                 lapply(1:2,function(s)fast_em(y,reduce_m22(m22,s)))))
 g3<-best_fit(list(g3,fast_em(y,expand_mixture_fit(m21),rep(1L,3))))
 m22<-best_fit(list(m22,fast_em(y,duplicate_singleton(m21))))
 # Refine until the declared stopping rule holds (not merely an iteration cap).
 finish<-function(f){for(i in 1:80){if(relative_gain(f)<=1e-9)break;f<-fast_em(y,f)};stopifnot(relative_gain(f)<=1.01e-9);f}
 fs<-lapply(list(G2=g2,G3=g3,M22=m22,M21=m21),finish)
 # The exact refinement is retained as a feasible G3 start.
 fs$G3<-finish(best_fit(list(fs$G3,fast_em(y,expand_mixture_fit(fs$M21),rep(1L,3)))))
 stopifnot(fs$G3$log_likelihood>=fs$M21$log_likelihood-1e-6,
           fs$M22$log_likelihood>=fs$M21$log_likelihood-1e-6)
 fs
}
