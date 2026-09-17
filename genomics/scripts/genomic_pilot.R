#!/usr/bin/env Rscript
library(Rcpp)
library(data.table)
library(parallel)
library(ggplot2)
root <- 'genomics'
base <- 'numerics'
source(file.path(base,'R/hmm_em.R'))
source(file.path(base,'R/mixture_hmm_em.R'))
source(file.path(base,'R/revision_validation.R'))
sourceCpp(file.path(base,'R/forward_backward.cpp'))
forward_backward_core <- forward_backward_compiled
out <- file.path(root,'results/genomic_pilot')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- fread(file.path(root,'data/genomics/development_chr1.csv'))
y <- d$beta_corrected
lengths <- as.integer(rle(d$sequence)$lengths)
stopifnot(sum(lengths)==length(y),all(is.finite(y)),all(d$chromosome==1L))
specifications <- list(G2=c(1,1),M21=c(2,1),M221=c(2,2,1))

refine <- function(f,counts,max_iter=2500,tol=1e-8) {
  mixture_em_loop(y,f$initial,f$transition,f$weights,f$means,f$sds,
                  counts,max_iter,tol,.01,lengths=lengths)
}
finalize <- function(f,counts) {
  f$K<-length(counts);f$M<-max(counts);f$component_counts<-as.integer(counts)
  f$n_parameters<-mixture_n_parameters(f)
  f$viterbi<-viterbi_mixture(y,f$initial,f$transition,f$weights,f$means,f$sds,lengths)
  f$BIC <- -2*f$log_likelihood+f$n_parameters*log(length(y))
  f$converged <- relative_gain(f)<=1e-8
  f
}
fit_model <- function(name,counts,nested=NULL) {
  set.seed(20260910+sum(counts)*100+length(counts))
  candidates <- lapply(seq_len(40),function(i)
    mixture_em_start(y,length(counts),max(counts),counts,100,1e-5,.01,lengths))
  # A variance-regime Gaussian model is a serious alternative to separated means.
  if(name=='G2')for(stay in c(.95,.995,.9995)) {
    p<-list(initial=c(.5,.5),transition=matrix(c(stay,1-stay,1-stay,stay),2),
            weights=matrix(1,2,1),means=matrix(.5,2,1),sds=matrix(c(.09,.32),2,1))
    candidates[[length(candidates)+1]]<-refine(p,counts,100,1e-5)
  }
  # Scientific starts pair low/high allele fractions within an imbalance state.
  # These supplement random starts; selection is by likelihood only.
  if(name %in% c('M21','M221')) for(stay in c(.95,.995,.9995)) {
    k<-length(counts);m<-max(counts)
    p<-list(initial=rep(1/k,k),transition=matrix((1-stay)/(k-1),k,k),
            weights=matrix(0,k,m),means=matrix(.5,k,m),sds=matrix(.08,k,m))
    diag(p$transition)<-stay
    for(s in seq_len(k))p$weights[s,seq_len(counts[s])]<-1/counts[s]
    p$means[1,]<-c(.1,.9)
    if(name=='M221')p$means[2,]<-c(.33,.67)
    candidates[[length(candidates)+1]]<-refine(p,counts,100,1e-5)
  }
  ll<-vapply(candidates,`[[`,numeric(1),'log_likelihood')
  selected<-head(order(ll,decreasing=TRUE),8)
  refined<-lapply(candidates[selected],refine,counts=counts)
  if(!is.null(nested)) {
    expanded<-expand_mixture_fit(nested)
    check<-forward_backward_mixture(y,expanded$initial,expanded$transition,
      expanded$weights,expanded$means,expanded$sds,lengths)$log_likelihood
    stopifnot(abs(check-nested$log_likelihood)<1e-6*(1+abs(check)))
    refined[[length(refined)+1]]<-refine(expanded,counts)
  }
  refined_ll<-vapply(refined,`[[`,numeric(1),'log_likelihood')
  f<-finalize(refined[[which.max(refined_ll)]],counts)
  f$audit<-rbind(data.frame(stage='screen',index=seq_along(ll),log_likelihood=ll),
                 data.frame(stage='refine',index=seq_along(refined_ll),log_likelihood=refined_ll))
  saveRDS(f,file.path(out,paste0(name,'.rds')))
  f
}
fits <- mclapply(names(specifications),function(n)fit_model(n,specifications[[n]]),mc.cores=3)
names(fits)<-names(specifications)
adopt_nested<-function(f,p,counts,name,stage) {
  check<-forward_backward_mixture(y,p$initial,p$transition,p$weights,p$means,p$sds,lengths)$log_likelihood
  for(jitter in c(0,.015,-.015)) {
    q<-p
    for(k in which(counts==2))q$means[k,]<-q$means[k,]+c(-jitter,jitter)
    candidate<-finalize(refine(q,counts),counts)
    row<-data.frame(stage=stage,index=jitter,log_likelihood=candidate$log_likelihood)
    audit<-rbind(f$audit,row)
    if(candidate$log_likelihood>f$log_likelihood)f<-candidate
    f$audit<-audit
  }
  stopifnot(f$log_likelihood>=check-1e-5)
  saveRDS(f,file.path(out,paste0(name,'.rds')));f
}
# G2 is contained in M21 by copying one emission component.
p<-fits$G2;p$weights<-matrix(c(.5,.5,1,0),2,byrow=TRUE)
p$means<-p$means[,c(1,1),drop=FALSE];p$sds<-p$sds[,c(1,1),drop=FALSE]
fits$M21<-adopt_nested(fits$M21,p,c(2,1),'M21','nested_G2')
fits$G3<-fit_model('G3',rep(1,3),nested=fits$M21)
# G3 is contained in M221 by copying the first two emission components.
p<-fits$G3;p$weights<-matrix(c(.5,.5,.5,.5,1,0),3,byrow=TRUE)
p$means<-p$means[,c(1,1),drop=FALSE];p$sds<-p$sds[,c(1,1),drop=FALSE]
fits$M221<-adopt_nested(fits$M221,p,c(2,2,1),'M221','nested_G3')
# M21 is contained in M221 by duplicating the first macrostate.
p<-fits$M21;idx<-c(1,1,2);w<-c(.5,.5,1)
p$initial<-p$initial[idx]*w
p$transition<-p$transition[idx,idx]*matrix(w,3,3,byrow=TRUE)
for(key in c('weights','means','sds'))p[[key]]<-p[[key]][idx,,drop=FALSE]
fits$M221<-adopt_nested(fits$M221,p,c(2,2,1),'M221','nested_M21')
fits$G5<-fit_model('G5',rep(1,5),nested=fits$M221)
stopifnot(fits$G3$log_likelihood>=fits$M21$log_likelihood-1e-5,
          fits$G5$log_likelihood>=fits$M221$log_likelihood-1e-5)
metrics<-list();params<-list()
for(n in names(fits)) {
  f<-fits[[n]]
  metrics[[n]]<-data.frame(model=n,states=f$K,components=sum(f$component_counts),
    parameters=f$n_parameters,log_likelihood=f$log_likelihood,BIC=f$BIC,
    runs=sum(vapply(split(f$viterbi,d$sequence),function(z)length(rle(z)$values),integer(1))),
    converged=f$converged,relative_gain=relative_gain(f),floor_hits=sum(f$sds<=.01*(1+1e-6)),
    minimum_transition=min(f$transition))
  for(k in seq_len(f$K))for(j in seq_len(f$component_counts[k]))params[[length(params)+1]]<-
    data.frame(model=n,state=k,component=j,weight=f$weights[k,j],mean=f$means[k,j],
               sd=f$sds[k,j],self_transition=f$transition[k,k],occupancy=mean(f$posterior[,k]))
  z<-copy(d);z[,state:=f$viterbi];fwrite(z,file.path(out,paste0(n,'_decoded.csv')))
}
fwrite(rbindlist(metrics),file.path(out,'metrics.csv'))
fwrite(rbindlist(params),file.path(out,'parameters.csv'))
print(rbindlist(metrics));print(rbindlist(params))
z<-rbindlist(lapply(names(fits),function(n) {
  x<-copy(d);x[,model:=n];x[,state:=factor(fits[[n]]$viterbi)];x
}))
p<-ggplot(z,aes(position/1e6,beta_corrected,color=state))+
  geom_point(size=.25,alpha=.55)+facet_wrap(~model,ncol=1)+
  labs(x='Chromosome 1 position (Mb, hg18)',y='Corrected tumor B-allele fraction',color='State')+
  theme_bw(base_size=11)+theme(legend.position='none')
ggsave(file.path(out,'development_profiles.png'),p,width=10,height=10,dpi=160)
