#!/usr/bin/env Rscript
library(Rcpp)
library(data.table)
library(parallel)
sourceCpp('genomics/scripts/gaussian_mixture_hmm.cpp')
source('genomics/scripts/symmetric_genomic_hmm.R')
root<-'genomics'
d<-fread(file.path(root,'data/genomics/development_chr1.csv'))
y<-d$beta_corrected;start<-c(TRUE,diff(d$sequence)!=0)
out<-file.path(root,'results/genomic_symmetry_benchmark')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
fits<-mclapply(2:3,function(k) {
  set.seed(20260910+k)
  candidates<-lapply(seq_len(30),function(i) {
    stay<-runif(1,.95,.9999)
    p<-list(initial=rep(1/k,k),transition=matrix((1-stay)/(k-1),k,k),
            mu=sort(runif(k,0,.48)),sd=runif(k,.025,.12))
    diag(p$transition)<-stay
    sym_fit(y,start,p,max_iter=120,tolerance=1e-5)
  })
  idx<-head(order(vapply(candidates,`[[`,numeric(1),'log_likelihood'),decreasing=TRUE),6)
  refined<-lapply(candidates[idx],function(p)sym_fit(y,start,p,max_iter=3000))
  f<-refined[[which.max(vapply(refined,`[[`,numeric(1),'log_likelihood'))]]
  f$audit<-data.table(log_likelihood=vapply(refined,`[[`,numeric(1),'log_likelihood'),
                     converged=vapply(refined,`[[`,logical(1),'converged'))
  saveRDS(f,file.path(out,paste0('F',k,'.rds')));f
},mc.cores=2)
names(fits)<-c('F2','F3')
rows<-rbindlist(lapply(names(fits),function(n) {
  f<-fits[[n]]
  data.table(model=n,states=length(f$mu),parameters=f$n_parameters,log_likelihood=f$log_likelihood,
    BIC=f$BIC,converged=f$converged,floor_hits=sum(f$sd<=.01*(1+1e-6)),
    runs=sum(vapply(split(f$viterbi,d$sequence),function(z)length(rle(z)$values),integer(1))))
}))
fwrite(rows,file.path(out,'metrics.csv'));print(rows)
parameters<-rbindlist(lapply(names(fits),function(n) {
  f<-fits[[n]];data.table(model=n,state=seq_along(f$mu),displacement=f$mu,sd=f$sd,
                         self_transition=diag(f$transition),occupancy=colMeans(f$posterior))
}))
fwrite(parameters,file.path(out,'parameters.csv'));print(parameters)
