#!/usr/bin/env Rscript
library(Rcpp)
library(data.table)
library(clue)
root<-'genomics'
sourceCpp(file.path(root,'scripts/gaussian_mixture_hmm.cpp'))
source(file.path(root,'scripts/symmetric_genomic_hmm.R'))
out<-file.path(root,'results/genomic_validation');dir.create(out,recursive=TRUE,showWarnings=FALSE)
development<-fread(file.path(root,'data/genomics/development_chr1.csv'))
datasets<-list(replication=fread(file.path(root,'data/genomics/replication_chr1.csv')),
               validation=fread(file.path(root,'data/genomics/validation_chr12_17.csv')))
confirmation_file<-file.path(root,'data/genomics/confirmation_chr18_22.csv')
if(file.exists(confirmation_file))datasets$confirmation<-fread(confirmation_file)
fits<-lapply(c('G2','M21','G3','G4','M221','G5'),function(n)readRDS(file.path(root,'results/genomic_pilot',paste0(n,'.rds'))))
names(fits)<-c('G2','M21','G3','G4','M221','G5')
for(n in c('F2','F3'))fits[[n]]<-readRDS(file.path(root,'results/genomic_symmetry_benchmark',paste0(n,'.rds')))
eval_fit<-function(d,f,symmetric) {
  start<-c(TRUE,diff(d$sequence)!=0)
  if(symmetric)return(sym_evaluate(d$beta_corrected,start,f,TRUE))
  e<-matrix(0,nrow(d),f$K)
  for(k in seq_len(f$K))for(j in seq_len(f$component_counts[k]))
    e[,k]<-e[,k]+f$weights[k,j]*dnorm(d$beta_corrected,f$means[k,j],f$sds[k,j])
  e<-pmax(e,1e-300)
  a<-fb_sequences_cpp(e,f$initial,f$transition,start)
  a$log_scores<-log(a$scales)
  a$viterbi<-viterbi_sequences_cpp(log(e),f$initial,f$transition,start)
  a
}
rows<-list();chromosomes<-list();replication<-list()
for(n in names(fits)) {
  f<-fits[[n]];sym<-startsWith(n,'F')
  baseline<-eval_fit(development,f,sym)
  stopifnot(abs(baseline$log_likelihood-f$log_likelihood)<1e-5*(1+abs(f$log_likelihood)))
  for(key in names(datasets)) {
    d<-copy(datasets[[key]]);e<-eval_fit(d,f,sym)
    d[,`:=`(state=e$viterbi,log_score=e$log_scores)]
    fwrite(d,file.path(out,paste0(key,'_',n,'.csv')))
    rows[[length(rows)+1]]<-data.table(dataset=key,model=n,n=nrow(d),log_score=sum(e$log_scores),
      mean_log_score=mean(e$log_scores),runs=sum(vapply(split(e$viterbi,d$sequence),function(z)length(rle(z)$values),integer(1))))
    x<-d[,.(n=.N,log_score=sum(log_score),mean_log_score=mean(log_score)),by=chromosome]
    model_name<-n
    x[,`:=`(dataset=key,model=model_name)];chromosomes[[length(chromosomes)+1]]<-x
    if(key=='replication') {
      reference<-development[,.(marker_index)];reference[,state_reference:=baseline$viterbi]
      common<-merge(d,reference,by='marker_index')
      replication[[n]]<-data.table(model=n,common_markers=nrow(common),state_agreement=mean(common$state==common$state_reference))
    }
  }
  # Fixed-parameter sensitivity to arbitrary allele designation.
  d<-copy(development);set.seed(20260910);flip<-sample(c(FALSE,TRUE),nrow(d),replace=TRUE)
  d[flip,beta_corrected:=1-beta_corrected]
  a<-eval_fit(d,f,sym)
  fwrite(data.table(model=n,log_likelihood_change=a$log_likelihood-baseline$log_likelihood,
    raw_state_agreement=mean(a$viterbi==baseline$viterbi)),file.path(out,paste0('allele_flip_',n,'.csv')))
}
fwrite(rbindlist(rows),file.path(out,'scores.csv'));print(rbindlist(rows))
fwrite(rbindlist(chromosomes),file.path(out,'chromosome_scores.csv'))
fwrite(rbindlist(replication),file.path(out,'technical_replication.csv'));print(rbindlist(replication))
# Emission-only correspondence between five ordinary states and M221 components.
m<-fits$M221;g<-fits$G5
idx<-do.call(rbind,lapply(seq_len(m$K),function(k)cbind(k,seq_len(m$component_counts[k]))))
means<-m$means[idx];sds<-m$sds[idx]
cost<-outer(seq_along(means),seq_len(g$K),Vectorize(function(i,j) {
  v<-sds[i]^2+g$sds[j,1]^2
  1-sqrt(2*sds[i]*g$sds[j,1]/v)*exp(-(means[i]-g$means[j,1])^2/(4*v))
}))
assignment<-as.integer(solve_LSAP(cost));mapping<-integer(g$K);mapping[assignment]<-idx[,1]
fwrite(data.table(mixture_state=idx[,1],mixture_component=idx[,2],ordinary_state=assignment,
                 hellinger_squared=cost[cbind(seq_along(means),assignment)]),file.path(out,'emission_alignment.csv'))
for(key in c('development',names(datasets))) {
  d<-if(key=='development')copy(development)else copy(datasets[[key]])
  a<-eval_fit(d,m,FALSE);b<-eval_fit(d,g,FALSE)
  collapsed<-mapping[b$viterbi]
  fwrite(data.table(dataset=key,collapsed_state_agreement=mean(collapsed==a$viterbi),
    mixture_runs=sum(vapply(split(a$viterbi,d$sequence),function(z)length(rle(z)$values),integer(1))),
    collapsed_ordinary_runs=sum(vapply(split(collapsed,d$sequence),function(z)length(rle(z)$values),integer(1)))),
    file.path(out,paste0('collapsed_',key,'.csv')))
  if(key=='development') {
    set.seed(20260910);flip<-sample(c(FALSE,TRUE),nrow(d),replace=TRUE)
    d[flip,beta_corrected:=1-beta_corrected]
    flipped<-eval_fit(d,g,FALSE)
    fwrite(data.table(model='G5_aggregated',state_agreement=mean(mapping[flipped$viterbi]==collapsed)),
           file.path(out,'allele_flip_G5_aggregated.csv'))
  }
}
