#!/usr/bin/env Rscript
library(testthat)
library(Rcpp)
library(dplyr)
testthat::test_dir("numerics/tests/testthat",reporter="summary",stop_on_failure=TRUE)
source("numerics/scripts/run_elk_application.R")
source("numerics/R/revision_validation.R")
Rcpp::sourceCpp("numerics/R/forward_backward.cpp")

test_that("compiled forward-backward matches enumeration and reference recursion", {
  set.seed(901)
  for(K in 1:4)for(T in c(1L,2L,4L)) {
    e<-matrix(runif(T*K,.001,2),T,K)
    initial<-normalize_probability(runif(K)); transition<-normalize_rows(matrix(runif(K*K),K))
    reference<-forward_backward_core(e,initial,transition)
    compiled<-forward_backward_compiled(e,initial,transition)
    expect_equal(compiled,reference,tolerance=1e-12)
    paths<-expand.grid(rep(list(seq_len(K)),T))
    probability<-apply(paths,1,function(path) {
      p<-initial[path[1]]*e[1,path[1]]
      if(T>1)for(t in 2:T)p<-p*transition[path[t-1],path[t]]*e[t,path[t]]
      p
    })
    expect_equal(exp(compiled$log_likelihood),sum(probability),tolerance=1e-12)
  }
})
forward_backward_core<-forward_backward_compiled
f<-readRDS("empirical/revision2/fits/elkelk-163_sd0.05_train158_M21.rds")
y<-read.csv("empirical/revision2/steps.csv") |> filter(ID=="elk-163") |> pull(log1p_step_km)
test_that("exact refinement preserves the whole observed likelihood", {
  expanded<-expand_mixture_fit(f)
  ll<-forward_backward(y,expanded$initial,expanded$transition,as.vector(expanded$means),as.vector(expanded$sds))$log_likelihood
  expect_equal(ll,f$log_likelihood,tolerance=1e-10)
  expect_equal(rowSums(expanded$transition),rep(1,3),tolerance=1e-12)
})
test_that("forecast scores telescope and future observations do not affect earlier forecasts", {
  train<-readRDS("empirical/revision2/fits/elkelk-163_sd0.05_train94_M21.rds")
  predicted<-sequential_predictions(y,train)
  changed<-y; changed[120:158]<-changed[120:158]+3
  expect_equal(predicted[1:119,],sequential_predictions(changed,train)[1:119,],tolerance=1e-12)
  ll<-function(z)forward_backward_mixture(z,train$initial,train$transition,train$weights,train$means,train$sds)$log_likelihood
  expect_equal(sum(predicted$log_score[95:126]),ll(y[1:126])-ll(y[1:94]),tolerance=1e-10)
})
test_that("all reported fits and validation values correspond to saved numerical objects", {
  summary<-read.csv("empirical/revision2/model_summary.csv")
  steps<-read.csv("empirical/revision2/steps.csv")
  pred<-read.csv("empirical/revision2/heldout_predictions.csv")
  for(i in seq_len(nrow(summary))) {
    r<-summary[i,]
    fit<-readRDS(sprintf("empirical/revision2/fits/elk%s_sd%s_train%s_%s.rds",r$ID,r$floor,r$n_train,r$model))
    series<-steps$log1p_step_km[steps$ID==r$ID]
    ll<-forward_backward_mixture(series[1:r$n_train],fit$initial,fit$transition,fit$weights,fit$means,fit$sds)$log_likelihood
    expect_equal(ll,r$log_likelihood,tolerance=1e-10)
    expect_equal(fit_bic(fit,r$n_train),r$BIC,tolerance=1e-10)
    expect_true(fit$converged)
    expect_true(all(diff(fit$log_likelihood_trace)>-1e-7))
    rows<-pred[pred$ID==r$ID&pred$model==r$model&pred$n_train==r$n_train,]
    if(nrow(rows)) {
      actual<-sequential_predictions(series,fit)[rows$step_index,]
      expect_equal(as.matrix(actual),as.matrix(rows[,c("log_score","pit","negative_mass")]),
                   tolerance=1e-10,ignore_attr=TRUE)
      expect_true(all(rows$step_index>r$n_train))
    }
  }
  for(ID in unique(summary$ID))for(model in unique(summary$model)) {
    rows<-summary[summary$ID==ID&summary$model==model&summary$fraction==1,]
    expect_true(all(diff(rows$log_likelihood[order(rows$floor)])<1e-4))
  }
  for(ID in unique(pred$ID))for(model in unique(pred$model)) {
    rows<-pred[pred$ID==ID&pred$model==model,]
    expect_equal(length(unique(rows$step_index)),nrow(rows))
  }
})
cat("All second-revision verification tests passed.\n")

source("numerics/scripts/test_optimization_audit.R")
