#!/usr/bin/env Rscript
library(Rcpp)
library(parallel)
library(dplyr)
Rcpp::sourceCpp('numerics/R/movement_likelihood.cpp')
source('numerics/R/movement_models.R')
out<-'empirical/application_review';data<-read.csv(file.path(out,'movement_data.csv'))
jobs<-expand.grid(ID=unique(data$ID),family=c('weibull','gamma'),fraction=c(.5,.6,.7,.8,.9),stringsAsFactors=FALSE)
models<-list(M21=c(2L,1L),H2=c(1L,1L),H3=rep(1L,3))
run<-function(j){
 job<-jobs[j,];whole<-data[data$ID==job$ID,];n<-floor(nrow(whole)*job$fraction);x<-head(whole,n)
 fits<-rows<-list()
 for(model in names(models)){
  name<-paste(job$ID,job$family,'joint',n,model,sep='_');path<-file.path(out,'validation_fits',paste0(name,'.rds'))
  previous<-readRDS(path)
  fit<-fit_movement(x,models[[model]],job$family,202609107L+100L*j+match(model,names(models)),
                    n_starts=80L,joint=TRUE,nested=if(model=='H3')fits$M21 else NULL,
                    extra_starts=list(previous$par))
  if(fit$nll<=previous$nll+1e-8){
    fit$previous_start_audit<-previous$audit
    fit$extra_start_improvement<-previous$nll-fit$nll
    saveRDS(fit,path)
  }else fit<-previous
  fits[[model]]<-fit
  rows[[model]]<-data.frame(ID=job$ID,family=job$family,n_train=n,model=model,
                          improvement=previous$nll-fit$nll,convergence=fit$convergence)
  cat(name,' improvement=',round(previous$nll-fit$nll,4),'\n')
 }
 bind_rows(rows)
}
results<-mclapply(seq_len(nrow(jobs)),run,mc.cores=5L,mc.preschedule=FALSE)
stopifnot(!any(vapply(results,inherits,logical(1),'try-error')))
write.csv(bind_rows(results),file.path(out,'validation_start_audit.csv'),row.names=FALSE)
