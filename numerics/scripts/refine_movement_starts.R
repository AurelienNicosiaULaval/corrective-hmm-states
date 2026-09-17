#!/usr/bin/env Rscript
library(Rcpp)
library(parallel)
library(dplyr)
Rcpp::sourceCpp('numerics/R/movement_likelihood.cpp')
source('numerics/R/movement_models.R')
out<-'empirical/application_review';data<-read.csv(file.path(out,'movement_data.csv'))
jobs<-expand.grid(ID=unique(data$ID),family=c('weibull','gamma'),stringsAsFactors=FALSE)
models<-list(M21=c(2L,1L),H2=c(1L,1L),H3=rep(1L,3))
run<-function(j){
 job<-jobs[j,];x<-data[data$ID==job$ID,];fits<-rows<-list()
 for(model in names(models)){
  name<-paste(job$ID,job$family,'joint',model,sep='_');path<-file.path(out,'fits',paste0(name,'.rds'))
  previous<-readRDS(path)
  candidates<-list.files(file.path(out,'sensitivity_fits'),pattern=paste0('^',job$ID,'_',job$family,'.*_',model,'[.]rds$'),full.names=TRUE)
  extra<-lapply(candidates,readRDS)
  extra<-lapply(extra[vapply(extra,function(f)length(f$par)==length(previous$par),logical(1))],`[[`,'par')
  fit<-fit_movement(x,models[[model]],job$family,202609105L+100L*j+match(model,names(models)),
                    n_starts=40L,joint=TRUE,nested=if(model=='H3')fits$M21 else NULL,
                    extra_starts=c(list(previous$par),extra))
  if(fit$nll<=previous$nll+1e-8){
    fit$previous_start_audit<-previous$audit
    fit$cross_start_improvement<-previous$nll-fit$nll
    saveRDS(fit,path)
  }else fit<-previous
  fits[[model]]<-fit
  rows[[model]]<-data.frame(ID=job$ID,family=job$family,model=model,improvement=previous$nll-fit$nll,
                          previous_BIC=previous$BIC,BIC=fit$BIC,convergence=fit$convergence)
 }
 bind_rows(rows)
}
results<-mclapply(seq_len(nrow(jobs)),run,mc.cores=4L,mc.preschedule=FALSE)
stopifnot(!any(vapply(results,inherits,logical(1),'try-error')))
write.csv(bind_rows(results),file.path(out,'cross_start_audit.csv'),row.names=FALSE)
