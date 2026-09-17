#!/usr/bin/env Rscript
library(Rcpp)
library(parallel)
library(dplyr)
source('numerics/scripts/run_elk_application.R')
source('numerics/R/revision_validation.R')
Rcpp::sourceCpp('numerics/R/forward_backward.cpp')
forward_backward_core<-forward_backward_compiled
out<-'empirical/application_review';dir.create(file.path(out,'gaussian_validation_fits'),showWarnings=FALSE)
data<-read.csv('empirical/revision2/steps.csv')
models<-list(M21=c(2L,1L),M22=c(2L,2L),G2=c(1L,1L),G3=rep(1L,3),G4=rep(1L,4))
jobs<-expand.grid(ID=unique(data$ID),fraction=c(.5,.6,.7,.8,.9),stringsAsFactors=FALSE)
run_job<-function(j){
 job<-jobs[j,];whole<-data[data$ID==job$ID,];n<-nrow(whole)
 n_train<-floor(job$fraction*n);end<-if(job$fraction==.9)n else floor((job$fraction+.1+1e-9)*n)
 y<-whole$log1p_step_km;fits<-rows<-list()
 for(m in seq_along(models)){
  model<-names(models)[m];key<-sprintf('elk%s_sd0.05_train%s_%s.rds',job$ID,n_train,model)
  old<-file.path('empirical/revision2/fits',key);path<-file.path(out,'gaussian_validation_fits',key)
  if(file.exists(old))fit<-readRDS(old)else if(file.exists(path))fit<-readRDS(path)else{
   nested<-if(model=='G3')fits$M21 else if(model=='G4')fits$M22 else NULL
   fit<-fit_audited(y[seq_len(n_train)],models[[m]],.05,202609103L+j*100L+m,nested)
   saveRDS(fit,path)
  }
  fits[[model]]<-fit
  p<-sequential_predictions(y,fit);use<-(n_train+1):end
  rows[[model]]<-cbind(whole[use,c('ID','step_index')],fraction=job$fraction,model=model,n_train=n_train,converged=fit$converged,p[use,])
  cat(key,' score=',round(sum(p$log_score[use]),2),'\n')
 }
 bind_rows(rows)
}
results<-mclapply(seq_len(nrow(jobs)),run_job,mc.cores=4L,mc.preschedule=FALSE)
stopifnot(!any(vapply(results,inherits,logical(1),'try-error')))
pred<-bind_rows(results)
write.csv(pred,file.path(out,'gaussian_validation_predictions.csv'),row.names=FALSE)
summary<-pred |> group_by(ID,model,fraction) |> summarise(n=n(),log_score=sum(log_score),converged=all(converged),.groups='drop')
write.csv(summary,file.path(out,'gaussian_validation_blocks.csv'),row.names=FALSE)
cat('Extended Gaussian validation complete.\n')
