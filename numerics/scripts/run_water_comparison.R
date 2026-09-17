#!/usr/bin/env Rscript
library(Rcpp)
library(parallel)
library(dplyr)
library(digest)
library(testthat)
Rcpp::sourceCpp('numerics/R/movement_likelihood.cpp')
source('numerics/R/movement_models.R')
out<-'empirical/application_review';dir.create(file.path(out,'water_fits'),showWarnings=FALSE)
data<-read.csv(file.path(out,'movement_data.csv'))
input<-tempfile(fileext='.txt')
download.file('https://ndownloader.figshare.com/files/5594180',input,quiet=TRUE)
stopifnot(digest(file=input,algo='sha256')=='36d423b2832aecd929f7d539da8779965b93a48fa6aa1321fdb7b194e68db2a5')
raw<-read.delim(input,check.names=FALSE);raw<-raw[raw$Individual!='',]
water<-bind_rows(lapply(split(raw,raw$Individual),function(x)data.frame(ID=x$Individual[1],
        interval=seq_len(nrow(x)-1),distance_water=x[['dist_water (meters)']][-nrow(x)]/1000)))
data<-left_join(data,water,by=c('ID','interval'))
stopifnot(!anyNA(data$distance_water))
write.csv(data,file.path(out,'movement_context.csv'),row.names=FALSE)

test_that('zero slopes reduce exactly to the homogeneous joint HMM',{
 x<-head(data[data$ID=='elk-287',],5);set.seed(44)
 p<-movement_initial(x,c(1L,1L),'weibull',TRUE)
 d<-movement_evaluate(p,x,c(1L,1L),'weibull',TRUE,TRUE)
 v<-movement_water_cpp(c(p,0,0),x$lower,x$upper,x$angle,1L,x$distance_water,TRUE)
 expect_equal(v$nll,d$nll,tolerance=1e-10)
 expect_equal(v$log_score,d$log_score,tolerance=1e-10)
 pp<-c(p,.7,-.3);v<-movement_water_cpp(pp,x$lower,x$upper,x$angle,1L,x$distance_water,TRUE)
 paths<-as.matrix(expand.grid(rep(list(1:2),5)))
 brute<-sum(apply(paths,1,function(s){
   z<-d$initial[s[1]]*exp(d$log_emission[1,s[1]])
   for(t in 2:5){q<-plogis(p[1+s[t-1]]+pp[length(p)+s[t-1]]*x$distance_water[t])
     z<-z*(if(s[t]==1)q else 1-q)*exp(d$log_emission[t,s[t]])}
   z
 }))
 expect_equal(-v$nll,log(brute),tolerance=1e-10)
 later<-x;later$distance_water[5]<-30
 changed<-movement_water_cpp(pp,x$lower,x$upper,x$angle,1L,later$distance_water,TRUE)
 expect_equal(v$log_score[1:4],changed$log_score[1:4],tolerance=1e-12)
})

jobs<-expand.grid(ID=unique(data$ID),family=c('weibull','gamma'),fraction=c(1,.5,.6,.7,.8,.9),stringsAsFactors=FALSE)
run<-function(j){
 job<-jobs[j,];whole<-data[data$ID==job$ID,];n<-nrow(whole);n_train<-floor(n*job$fraction)
 x<-head(whole,n_train);covariate<-log1p(x$distance_water)
 centre<-mean(covariate);spread<-sd(covariate);z<-(covariate-centre)/spread
 basepath<-if(job$fraction==1)file.path(out,'fits',paste0(job$ID,'_',job$family,'_joint_H2.rds'))else
   file.path(out,'validation_fits',paste0(job$ID,'_',job$family,'_joint_',n_train,'_H2.rds'))
 base<-readRDS(basepath)
 path<-file.path(out,'water_fits',paste0(job$ID,'_',job$family,'_',n_train,'.rds'))
 objective<-function(p)movement_water_cpp(p,x$lower,x$upper,x$angle,movement_families[[job$family]],z,FALSE)$nll
 bounds<-movement_bounds(c(1L,1L),job$family,TRUE)
 lower<-c(bounds$lower,-10,-10);upper<-c(bounds$upper,10,10)
 if(file.exists(path))fit<-readRDS(path)else{
   set.seed(202609104L+j)
   starts<-c(list(c(base$par,0,0)),lapply(1:29,function(k)c(movement_initial(x,c(1L,1L),job$family,TRUE),rnorm(2))))
   optimize<-function(p,iterations,tol)nlminb(p,objective,lower=lower,upper=upper,
                         control=list(iter.max=iterations,eval.max=iterations*4,rel.tol=tol))
   screen<-lapply(starts,optimize,iterations=100,tol=1e-4)
   order<-order(vapply(screen,`[[`,numeric(1),'objective'))
   refined<-lapply(screen[head(order,8)],function(s)optimize(s$par,3000,1e-10))
   best<-refined[[which.min(vapply(refined,`[[`,numeric(1),'objective'))]]
   if(best$convergence!=0){
     alternative<-optim(best$par,objective,method='L-BFGS-B',lower=lower,upper=upper,
                        control=list(maxit=5000,factr=1e6,ndeps=rep(1e-5,length(best$par))))
     restart<-optimize(alternative$par,5000,1e-10)
     if(restart$objective<=best$objective+1e-7 && restart$convergence==0)best<-restart
   }
   fit<-list(par=best$par,nll=best$objective,convergence=best$convergence,message=best$message,
             family=job$family,n_intervals=n_train,n_parameters=length(best$par),
             centre=centre,spread=spread,screen=screen,refined=refined,seed=202609104L+j,
             BIC=2*best$objective+length(best$par)*log(sum(is.finite(x$rate)|is.finite(x$angle))))
   saveRDS(fit,path)
 }
 # Restore nesting after any refinement of the training-only homogeneous fit.
 if(fit$nll>base$nll+1e-8){
   z<-nlminb(c(base$par,0,0),objective,lower=lower,upper=upper,
             control=list(iter.max=5000,eval.max=20000,rel.tol=1e-9))
   if(z$objective<fit$nll){
     fit$par<-z$par;fit$nll<-z$objective;fit$convergence<-z$convergence;fit$message<-z$message
     fit$BIC<-2*fit$nll+fit$n_parameters*log(sum(is.finite(x$rate)|is.finite(x$angle)))
     fit$nesting_refresh<-z;saveRDS(fit,path)
   }
 }
 details<-movement_water_cpp(fit$par,whole$lower,whole$upper,whole$angle,movement_families[[job$family]],
                         (log1p(whole$distance_water)-fit$centre)/fit$spread,TRUE)
 use<-if(job$fraction==1)seq_len(n)else (n_train+1):(if(job$fraction==.9)n else floor((job$fraction+.1+1e-9)*n))
 cat(job$ID,job$family,n_train,'BIC=',round(fit$BIC,2),' convergence=',fit$convergence,'\n')
 cbind(job,data.frame(n_train=n_train,BIC=fit$BIC,n_parameters=fit$n_parameters,convergence=fit$convergence,
                     score=sum(details$log_score[use]),homogeneous_BIC=base$BIC))
}
results<-mclapply(seq_len(nrow(jobs)),run,mc.cores=3L,mc.preschedule=FALSE)
stopifnot(!any(vapply(results,inherits,logical(1),'try-error')))
write.csv(bind_rows(results),file.path(out,'water_comparison.csv'),row.names=FALSE)
