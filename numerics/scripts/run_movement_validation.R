#!/usr/bin/env Rscript
library(Rcpp)
library(parallel)
library(dplyr)
Rcpp::sourceCpp('numerics/R/movement_likelihood.cpp')
source('numerics/R/movement_models.R')
out<-'empirical/application_review';dir.create(file.path(out,'validation_fits'),showWarnings=FALSE)
data<-read.csv(file.path(out,'movement_data.csv'))
models<-list(M21=c(2L,1L),H2=c(1L,1L),H3=rep(1L,3))
joint<-'joint'%in%commandArgs(trailingOnly=TRUE)
jobs<-expand.grid(ID=unique(data$ID),family=names(movement_families),fraction=c(.5,.6,.7,.8,.9),stringsAsFactors=FALSE)
run_job<-function(j){
 job<-jobs[j,];whole<-data[data$ID==job$ID,];n<-nrow(whole)
 n_train<-floor(job$fraction*n);end<-if(job$fraction==.9)n else floor((job$fraction+.1+1e-9)*n)
 train<-whole[seq_len(n_train),];fits<-list();rows<-list()
 for(m in seq_along(models)){
  model<-names(models)[m];name<-paste(job$ID,job$family,if(joint)'joint'else'rate',n_train,model,sep='_')
  path<-file.path(out,'validation_fits',paste0(name,'.rds'))
  if(file.exists(path))fit<-readRDS(path)else{
   nested<-if(model=='H3')fits$M21 else NULL
   fit<-fit_movement(train,models[[m]],job$family,202609102L+j*100L+m,n_starts=20L,joint=joint,nested=nested)
   saveRDS(fit,path)
  }
  fits[[model]]<-fit
  pred<-movement_predictions(fit,whole)
  use<-(n_train+1):end
  rows[[model]]<-cbind(whole[use,c('ID','interval','rate','angle')],family=job$family,joint=joint,
    fraction=job$fraction,model=model,n_train=n_train,convergence=fit$convergence,pred[use,])
  cat(name,' score=',round(sum(pred$log_score[use]),2),' convergence=',fit$convergence,'\n')
 }
 bind_rows(rows)
}
results<-mclapply(seq_len(nrow(jobs)),run_job,mc.cores=3L,mc.preschedule=FALSE)
stopifnot(!any(vapply(results,inherits,logical(1),'try-error')))
pred<-bind_rows(results)
write.csv(pred,file.path(out,if(joint)'joint_validation_predictions.csv'else'rate_validation_predictions.csv'),row.names=FALSE)
summary<-pred |> group_by(ID,family,model,joint,fraction) |> summarise(n=n(),n_rates=sum(is.finite(rate)),
  log_score=sum(log_score),rate_log_score=sum(rate_log_score,na.rm=TRUE),converged=all(convergence==0),.groups='drop')
write.csv(summary,file.path(out,if(joint)'joint_validation_blocks.csv'else'rate_validation_blocks.csv'),row.names=FALSE)
cat('Temporal validation complete.\n')
