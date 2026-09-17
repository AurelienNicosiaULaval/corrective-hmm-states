#!/usr/bin/env Rscript
library(Rcpp)
library(testthat)
library(dplyr)
Rcpp::sourceCpp('numerics/R/movement_likelihood.cpp')
source('numerics/R/movement_models.R')
out<-'empirical/application_review';data<-read.csv(file.path(out,'movement_context.csv'))
test_that('all saved primary fits match their actual observations and BIC',{
 for(folder in c('fits','validation_fits'))for(path in list.files(file.path(out,folder),pattern='[.]rds$',full.names=TRUE)){
  f<-readRDS(path);id<-sub('_.*','',basename(path));x<-head(data[data$ID==id,],f$n_intervals)
  d<-movement_evaluate(f$par,x,f$counts,f$family,f$joint,TRUE,f$shared_angle)
  expect_equal(d$nll,f$nll,tolerance=1e-8,info=path)
  expect_equal(f$BIC,2*f$nll+f$n_parameters*log(sum(is.finite(x$rate)|(f$joint&is.finite(x$angle)))),tolerance=1e-8)
  expect_equal(rowSums(d$transition),rep(1,length(f$counts)),tolerance=1e-12)
  expect_equal(length(f$par),f$n_parameters)
 }
})
test_that('every saved forecast is obtained from its own training-only fit',{
 for(type in c('joint','rate')){
  pred<-read.csv(file.path(out,paste0(type,'_validation_predictions.csv')))
  groups<-split(pred,interaction(pred$ID,pred$family,pred$n_train,pred$model,drop=TRUE))
  for(z in groups){
   path<-file.path(out,'validation_fits',paste0(z$ID[1],'_',z$family[1],'_',type,'_',z$n_train[1],'_',z$model[1],'.rds'))
   f<-readRDS(path);x<-data[data$ID==z$ID[1],]
   expected<-movement_predictions(f,x)[z$interval,]
   expect_true(all(z$interval>f$n_intervals))
   expect_equal(z$log_score,expected$log_score,tolerance=1e-8,info=path)
   expect_equal(z$cdf_lower,expected$cdf_lower,tolerance=1e-8)
   expect_equal(z$cdf_upper,expected$cdf_upper,tolerance=1e-8)
  }
  distinct_keys<-pred |> distinct(ID,family,model,interval)
  expect_equal(nrow(distinct_keys),nrow(pred))
 }
})
test_that('ordinary models satisfy their nesting likelihood inequalities',{
 for(folder in c('fits','validation_fits')){
  files<-list.files(file.path(out,folder),pattern='_M21[.]rds$',full.names=TRUE)
  for(path in files){
   mix<-readRDS(path);h<-readRDS(sub('_M21[.]rds$','_H3.rds',path))
   expect_true(h$nll<=mix$nll+1e-5,info=path)
  }
 }
})
test_that('water fits use the same observations and contain their homogeneous alternatives',{
 tab<-read.csv(file.path(out,'water_comparison.csv'))
 for(j in seq_len(nrow(tab))){
  row<-tab[j,];x<-head(data[data$ID==row$ID,],row$n_train)
  path<-file.path(out,'water_fits',paste0(row$ID,'_',row$family,'_',row$n_train,'.rds'))
  f<-readRDS(path)
  d<-movement_water_cpp(f$par,x$lower,x$upper,x$angle,movement_families[[f$family]],
                        (log1p(x$distance_water)-f$centre)/f$spread,TRUE)
  expect_equal(d$nll,f$nll,tolerance=1e-8)
  expect_equal(f$centre,mean(log1p(x$distance_water)),tolerance=1e-12)
  expect_equal(f$spread,sd(log1p(x$distance_water)),tolerance=1e-12)
  baseline<-if(row$fraction==1)file.path(out,'fits',paste0(row$ID,'_',row$family,'_joint_H2.rds'))else
    file.path(out,'validation_fits',paste0(row$ID,'_',row$family,'_joint_',row$n_train,'_H2.rds'))
  expect_true(f$nll<=readRDS(baseline)$nll+1e-5,info=path)
 }
})
test_that('resolution and angular sensitivity fits match the stated observations',{
 tab<-read.csv(file.path(out,'sensitivity_summary.csv'))
 for(j in seq_len(nrow(tab))){
  r<-tab[j,];x<-data[data$ID==r$ID,]
  if(grepl('coarse',r$scenario)){
    threshold<-if(r$scenario=='coarse_005').05 else .1
    small<-is.finite(x$rate)&x$rate<threshold-1e-8
    x$lower[small]<-0;x$upper[small]<-threshold-.005
    x$angle[small|c(TRUE,head(small,-1))]<-NA_real_
  }
  path<-file.path(out,'sensitivity_fits',paste0(r$ID,'_',r$family,'_',r$scenario,'_',r$model,'.rds'))
  f<-readRDS(path)
  z<-movement_evaluate(f$par,x,f$counts,f$family,TRUE,FALSE,f$shared_angle)$nll
  expect_equal(z,f$nll,tolerance=1e-8,info=path)
  expect_equal(r$BIC,2*z+f$n_parameters*log(sum(is.finite(x$rate)|is.finite(x$angle))),tolerance=1e-8)
 }
})
cat('Application review results verification passed.\n')
