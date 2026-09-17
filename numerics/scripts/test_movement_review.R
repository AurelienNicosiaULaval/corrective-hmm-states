#!/usr/bin/env Rscript
library(Rcpp)
library(testthat)
Rcpp::sourceCpp('numerics/R/movement_likelihood.cpp')
source('numerics/R/movement_models.R')

test_that('one-state interval probabilities match independent R distribution calculations', {
 d<-data.frame(lower=c(0,.005,.05,.5,2,NA),upper=c(.005,.015,.06,.51,2.01,NA),angle=NA_real_,rate=c(0,.01,.055,.505,2.005,NA))
 for(family in names(movement_families)){
   par<-if(family=='truncated_normal')c(.3,log(.4))else c(log(.7),log(1.2))
   got<-movement_evaluate(par,d,1L,family,detail=TRUE)
   cdf<-switch(family,gamma=function(x)pgamma(x,shape=exp(par[2]),scale=exp(par[1]-par[2])),
     weibull=function(x)pweibull(x,shape=exp(par[2]),scale=exp(par[1])),
     lognormal=function(x)plnorm(x,meanlog=par[1],sdlog=exp(par[2])),
     truncated_normal=function(x)(pnorm(log1p(x),par[1],exp(par[2]))-pnorm(0,par[1],exp(par[2])))/pnorm(0,par[1],exp(par[2]),lower.tail=FALSE))
   expected<-log(cdf(d$upper)-cdf(d$lower));expected[6]<-0
   expect_equal(as.numeric(got$log_emission),expected,tolerance=1e-10)
   expect_equal(got$nll,-sum(expected),tolerance=1e-10)
 }
})
test_that('forward likelihood equals explicit path enumeration and exact component refinement', {
 set.seed(35)
 d<-data.frame(lower=c(.05,.5,1),upper=c(.06,.51,1.01),angle=c(.3,-2,.2),rate=c(.055,.505,1.005))
 for(joint in c(FALSE,TRUE)){
  counts<-c(2L,1L);par<-movement_initial(d,counts,'gamma',joint)
  fit<-list(par=par,counts=counts,family='gamma',joint=joint,shared_angle=FALSE)
  a<-movement_evaluate(par,d,counts,'gamma',joint,TRUE)
  paths<-expand.grid(rep(list(1:2),3))
  prob<-apply(paths,1,function(z)a$initial[z[1]]*exp(sum(a$log_emission[cbind(1:3,z)]))*a$transition[z[1],z[2]]*a$transition[z[2],z[3]])
  expect_equal(-a$nll,log(sum(prob)),tolerance=1e-12)
  expanded<-movement_expansion(fit,d)
  b<-movement_evaluate(expanded,d,rep(1L,3),'gamma',joint,TRUE)
  expect_equal(a$nll,b$nll,tolerance=1e-12)
  expect_equal(rowSums(a$transition),rep(1,2),tolerance=1e-12)
 }
})
test_that('scores use past observations only and zero is treated as an interval', {
 d<-data.frame(lower=c(0,.005,.2,.5),upper=c(.005,.015,.21,.51),angle=c(NA,.3,-2,0),rate=c(0,.01,.205,.505))
 par<-movement_initial(d,c(2L,1L),'lognormal',TRUE)
 a<-movement_evaluate(par,d,c(2L,1L),'lognormal',TRUE,TRUE)
 changed<-d;changed$lower[4]<-5;changed$upper[4]<-5.01
 b<-movement_evaluate(par,changed,c(2L,1L),'lognormal',TRUE,TRUE)
 expect_equal(a$log_score[1:3],b$log_score[1:3],tolerance=1e-12)
 expect_true(all(is.finite(a$log_score)))
 expect_true(all(exp(a$log_components[1,])>0))
})
test_that('binned forecast CDFs agree with rate probabilities and scores telescope', {
 d<-data.frame(lower=c(0,.005,.2,.5),upper=c(.005,.015,.21,.51),angle=c(NA,.3,-2,0),rate=c(0,.01,.205,.505))
 for(family in names(movement_families)){
  fit<-list(par=movement_initial(d,c(2L,1L),family,TRUE),counts=c(2L,1L),family=family,joint=TRUE,shared_angle=FALSE)
  p<-movement_predictions(fit,d)
  expect_equal(exp(p$rate_log_score),p$cdf_upper-p$cdf_lower,tolerance=1e-10)
  a<-movement_evaluate(fit$par,d,fit$counts,family,TRUE)$nll
  b<-movement_evaluate(fit$par,d[1:2,],fit$counts,family,TRUE)$nll
  expect_equal(sum(p$log_score[3:4]),b-a,tolerance=1e-10)
 }
})
cat('Movement likelihood verification passed.\n')
