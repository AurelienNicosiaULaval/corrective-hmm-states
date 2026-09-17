#!/usr/bin/env Rscript
library(Rcpp)
library(data.table)
library(parallel)
root<-normalizePath('genomics');.libPaths(c(file.path(root,'library'),.libPaths()))
library(aroma.light)
base<-'numerics'
source(file.path(base,'R/hmm_em.R'));source(file.path(base,'R/mixture_hmm_em.R'))
source(file.path(base,'R/revision_validation.R'));sourceCpp(file.path(base,'R/forward_backward.cpp'))
forward_backward_core<-forward_backward_compiled
p<-file.path(root,'data/genomics');out<-file.path(root,'results/genomic_sensitivity')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
reference<-fread(file.path(root,'results/genomic_pilot/M221_decoded.csv'))
baseline<-readRDS(file.path(root,'results/genomic_pilot/M221.rds'))
coords<-as.data.table(readRDS(file.path(p,'snp6_coordinates_hg18.rds')))
t<-readRDS(file.path(p,'GSM337641_signals.rds'));n<-readRDS(file.path(p,'GSM337662_signals.rds'))
coords[,`:=`(marker_index=.I,beta_tumor=t$fracB,beta_normal=n$fracB)]
coords<-coords[chromosome==1 & is.finite(position) & position>0 & is.finite(beta_tumor) & is.finite(beta_normal)]
settings<-data.table(name=c('floor_005','floor_020','thin_2','thin_5','normal_025','normal_035'),
                    floor=c(.005,.02,.01,.01,.01,.01),thin=c(1,1,2,5,1,1),
                    low=c(.3,.3,.3,.3,.25,.35))
results<-mclapply(seq_len(nrow(settings)),function(i) {
  setting<-settings[i];lower<-setting$low
  d<-copy(coords[beta_normal>=lower & beta_normal<=1-lower]);setorder(d,position,marker_index)
  d<-d[!duplicated(position)];d[,sequence:=cumsum(c(TRUE,diff(position)>1e6))]
  d[,within_sequence:=seq_len(.N),by=sequence]
  d<-d[(within_sequence-1)%%setting$thin==0]
  d[,sequence:=cumsum(c(TRUE,diff(position)>1e6))]
  d[,beta_corrected:=as.numeric(normalizeTumorBoost(beta_tumor,beta_normal,
     muN=rep(.5,.N),preserveScale=FALSE,flavor='v4'))]
  y<-d$beta_corrected;lengths<-as.integer(rle(d$sequence)$lengths);counts<-c(2L,2L,1L)
  floor<-setting$floor;set.seed(20260910+i)
  refine<-function(f,it=2000,tol=1e-8)mixture_em_loop(y,f$initial,f$transition,
     f$weights,f$means,f$sds,counts,it,tol,floor,lengths)
  initial<-baseline
  if(setting$thin>1)initial$transition<-Reduce(`%*%`,rep(list(baseline$transition),setting$thin))
  candidates<-lapply(seq_len(10),function(i)mixture_em_start(y,3,2,counts,100,1e-5,floor,lengths))
  top<-head(order(vapply(candidates,`[[`,numeric(1),'log_likelihood'),decreasing=TRUE),3)
  fits<-c(lapply(candidates[top],refine),list(refine(initial)))
  f<-fits[[which.max(vapply(fits,`[[`,numeric(1),'log_likelihood'))]]
  f$K<-3;f$M<-2;f$component_counts<-counts
  f$viterbi<-viterbi_mixture(y,f$initial,f$transition,f$weights,f$means,f$sds,lengths)
  # Map paired states by their separation; keep the one-component state separate.
  order<-order(abs(f$means[1:2,1]-f$means[1:2,2]),decreasing=TRUE)
  mapping<-integer(3);mapping[order]<-1:2;mapping[3]<-3
  d[,state:=mapping[f$viterbi]]
  common<-merge(d[,.(marker_index,state)],reference[,.(marker_index,state_reference=state)],by='marker_index')
  saveRDS(f,file.path(out,paste0(setting$name,'.rds')));fwrite(d,file.path(out,paste0(setting$name,'_decoded.csv')))
  data.table(setting=setting$name,n=nrow(d),floor=floor,common_markers=nrow(common),
    state_agreement=mean(common$state==common$state_reference),log_likelihood=f$log_likelihood,
    converged=relative_gain(f)<=1e-8,floor_hits=sum(f$sds<=floor*(1+1e-6)),
    runs=sum(vapply(split(d$state,d$sequence),function(z)length(rle(z)$values),integer(1))))
},mc.cores=4)
fwrite(rbindlist(results),file.path(out,'summary.csv'));print(rbindlist(results))
