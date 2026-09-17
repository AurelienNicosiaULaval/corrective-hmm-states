#!/usr/bin/env Rscript
library(Rcpp)
library(parallel)
library(dplyr)
library(digest)
library(jsonlite)
Rcpp::sourceCpp('numerics/R/movement_likelihood.cpp')
source('numerics/R/movement_models.R')
out<-'empirical/application_review';dir.create(file.path(out,'fits'),recursive=TRUE,showWarnings=FALSE)
input<-tempfile(fileext='.txt')
download.file('https://ndownloader.figshare.com/files/5594180',input,quiet=TRUE)
stopifnot(digest(file=input,algo='sha256')=='36d423b2832aecd929f7d539da8779965b93a48fa6aa1321fdb7b194e68db2a5')
raw<-read.delim(input,check.names=FALSE);raw<-raw[raw$Individual!='',]
data<-prepare_movement_data(raw);write.csv(data,file.path(out,'movement_data.csv'),row.names=FALSE)
models<-list(M21=c(2L,1L),M22=c(2L,2L),H2=c(1L,1L),H3=rep(1L,3),H4=rep(1L,4))
arguments<-commandArgs(trailingOnly=TRUE)
joint<-'joint'%in%arguments
jobs<-expand.grid(ID=unique(data$ID),family=names(movement_families),stringsAsFactors=FALSE)
run_job<-function(j){
  job<-jobs[j,];d<-data[data$ID==job$ID,];fits<-list()
  for(m in seq_along(models)){
    model<-names(models)[m];name<-paste(job$ID,job$family,if(joint)'joint'else'rate',model,sep='_')
    path<-file.path(out,'fits',paste0(name,'.rds'))
    if(file.exists(path))fit<-readRDS(path)else{
      nested<-if(model=='H3')fits$M21 else if(model=='H4')fits$M22 else NULL
      fit<-fit_movement(d,models[[m]],job$family,202609101L+j*100L+m,n_starts=40L,joint=joint,nested=nested)
      saveRDS(fit,path)
    }
    fits[[model]]<-fit
    cat(name,' BIC=',round(fit$BIC,2),' convergence=',fit$convergence,'\n')
  }
  fits
}
results<-mclapply(seq_len(nrow(jobs)),run_job,mc.cores=4L,mc.preschedule=FALSE)
stopifnot(!any(vapply(results,inherits,logical(1),'try-error')))
rows<-list()
for(j in seq_len(nrow(jobs)))for(model in names(models)){
 f<-results[[j]][[model]]
 rows[[length(rows)+1]]<-cbind(jobs[j,],data.frame(model=model,joint=joint,n=f$n_intervals,
   n_parameters=f$n_parameters,log_likelihood=-f$nll,BIC=f$BIC,convergence=f$convergence,
   parameters_at_bounds=length(f$boundary_indices)))
}
write.csv(bind_rows(rows),file.path(out,if(joint)'joint_summary.csv'else'rate_summary.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
