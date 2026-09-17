#!/usr/bin/env Rscript
library(Rcpp)
library(dplyr)
library(parallel)
Rcpp::sourceCpp('numerics/R/movement_likelihood.cpp')
Rcpp::sourceCpp('numerics/R/forward_backward.cpp')
source('numerics/R/movement_models.R')
out <- 'empirical/application_review'
dir.create(file.path(out,'sensitivity_fits'),showWarnings=FALSE)
data <- read.csv(file.path(out,'movement_data.csv'))
models <- list(M21=c(2L,1L),H2=c(1L,1L),H3=rep(1L,3))
jobs <- expand.grid(ID=unique(data$ID),family=c('weibull','gamma'),
                    scenario=c('coarse_005','coarse_010','shared_angle','rho_090','rho_098'),stringsAsFactors=FALSE)
run <- function(j) {
 job<-jobs[j,];x<-data[data$ID==job$ID,];base_x<-x
 if(grepl('coarse',job$scenario)){
   threshold<-if(job$scenario=='coarse_005').05 else .1
   small<-is.finite(x$rate)&x$rate<threshold-1e-8
   x$lower[small]<-0;x$upper[small]<-threshold-.005
   # Angles involving either of the two smallest steps are uninformative at this resolution.
   imprecise<-small | c(TRUE,head(small,-1))
   x$angle[imprecise]<-NA_real_
 }
 shared<-job$scenario=='shared_angle'
 rho_max<-if(job$scenario=='rho_090').9 else if(job$scenario=='rho_098').98 else .95
 fits<-list();rows<-list()
 for(model in names(models)){
   name<-paste(job$ID,job$family,job$scenario,model,sep='_');path<-file.path(out,'sensitivity_fits',paste0(name,'.rds'))
   baseline<-readRDS(file.path(out,'fits',paste0(job$ID,'_',job$family,'_joint_',model,'.rds')))
   if(file.exists(path))fit<-readRDS(path)else{
     extra<-if(!shared||model!='M21')list(baseline$par)else list()
     # With shared angular laws the ordinary model still has one law per state.
     fit<-fit_movement(x,models[[model]],job$family,202609103L+j*100L+match(model,names(models)),
                       n_starts=30L,joint=TRUE,shared_angle=shared,nested=if(model=='H3')fits$M21 else NULL,
                       rho_max=rho_max,extra_starts=extra)
     saveRDS(fit,path)
   }
   # Refresh cached solutions from the final baseline and exact nested fit.
   bounds<-movement_bounds(models[[model]],job$family,TRUE,shared,rho_max)
   objective<-function(p)movement_evaluate(p,x,models[[model]],job$family,TRUE,FALSE,shared)$nll
   starts<-list(fit$par)
   if(length(baseline$par)==length(fit$par))starts<-c(starts,list(baseline$par))
   if(model=='H3')starts<-c(starts,list(movement_expansion(fits$M21,x)))
   for(start in starts){
     z<-nlminb(pmax(bounds$lower,pmin(bounds$upper,start)),objective,lower=bounds$lower,upper=bounds$upper,
              control=list(iter.max=2000,eval.max=6000,rel.tol=1e-9))
     if(z$objective<fit$nll){
       fit$par<-z$par;fit$nll<-z$objective;fit$convergence<-z$convergence;fit$message<-z$message
       fit$BIC<-2*fit$nll+fit$n_parameters*log(sum(is.finite(x$rate)|is.finite(x$angle)))
       fit$boundary_indices<-which(abs(fit$par-bounds$lower)<1e-4|abs(fit$par-bounds$upper)<1e-4)
       fit$baseline_refresh<-z
     }
   }
   saveRDS(fit,path)
   fits[[model]]<-fit
   if(model=='M21'){
     now<-movement_posteriors(fit,x);old<-movement_posteriors(baseline,base_x)
     agreement<-max(mean(now$viterbi==old$viterbi),mean(3-now$viterbi==old$viterbi))
     runs<-sum(diff(now$viterbi)!=0)+1L
   }else{agreement<-NA;runs<-NA}
   rows[[model]]<-cbind(job,data.frame(model=model,BIC=fit$BIC,log_likelihood=-fit$nll,
                       n_parameters=fit$n_parameters,convergence=fit$convergence,
                       agreement=agreement,decoded_runs=runs,bounds=length(fit$boundary_indices)))
   cat(name,round(fit$BIC,2),' convergence=',fit$convergence,'\n')
 }
 bind_rows(rows)
}
results<-mclapply(seq_len(nrow(jobs)),run,mc.cores=4L,mc.preschedule=FALSE)
stopifnot(!any(vapply(results,inherits,logical(1),'try-error')))
write.csv(bind_rows(results),file.path(out,'sensitivity_summary.csv'),row.names=FALSE)
