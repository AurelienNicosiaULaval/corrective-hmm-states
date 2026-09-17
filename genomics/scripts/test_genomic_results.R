#!/usr/bin/env Rscript
library(Rcpp)
library(data.table)
library(testthat)
source('numerics/R/hmm_em.R')
source('numerics/R/mixture_hmm_em.R')
source('numerics/R/revision_validation.R')
sourceCpp('numerics/R/forward_backward.cpp')
forward_backward_core <- forward_backward_compiled
sourceCpp('genomics/scripts/gaussian_mixture_hmm.cpp')
source('genomics/scripts/symmetric_genomic_hmm.R')
p <- 'genomics'
readfit <- function(n)readRDS(file.path(p,'results',if(startsWith(n,'F'))'genomic_symmetry_benchmark' else 'genomic_pilot',paste0(n,'.rds')))
d <- fread(file.path(p,'data/genomics/development_chr1.csv'))
start<-c(TRUE,diff(d$sequence)!=0);lengths<-as.integer(rle(d$sequence)$lengths)
metrics<-rbindlist(list(fread(file.path(p,'results/genomic_pilot/metrics.csv')),fread(file.path(p,'results/genomic_symmetry_benchmark/metrics.csv'))),fill=TRUE)
test_that('partition, source rows and saved fit summaries agree',{
 expect_equal(nrow(d),20697L);expect_equal(lengths,c(10322L,10375L))
 expect_true(all(d$beta_normal>=.3 & d$beta_normal<=.7))
 expect_equal(anyDuplicated(d[,.(chromosome,position)]),0L)
 expect_equal(unique(d$chromosome),1L)
 val<-fread(file.path(p,'data/genomics/validation_chr12_17.csv'))
 conf<-fread(file.path(p,'data/genomics/confirmation_chr18_22.csv'))
 expect_equal(unique(val$chromosome),12:17);expect_equal(unique(conf$chromosome),18:22)
 expect_length(intersect(d$marker_index,val$marker_index),0)
 expect_length(intersect(val$marker_index,conf$marker_index),0)
 for(n in c('G2','M21','G3','M221','G5')) {
  f<-readfit(n);expect_true(isTRUE(f$converged));expect_true(all(f$sds>=.01))
  expect_true(all(diff(f$log_likelihood_trace)>= -1e-6*(1+abs(f$log_likelihood))))
  expect_equal(f$BIC,-2*f$log_likelihood+f$n_parameters*log(nrow(d)),tolerance=1e-10)
  expect_equal(f$log_likelihood,metrics[model==n,log_likelihood],tolerance=1e-10)
  expect_equal(sum(f$initial),1,tolerance=1e-10);expect_equal(rowSums(f$transition),rep(1,f$K),tolerance=1e-10)
 }
})
test_that('component refinement reproduces multiple-sequence likelihood and inclusion',{
 for(pair in list(c('M21','G3'),c('M221','G5'))) {
  f<-readfit(pair[1]);g<-readfit(pair[2]);expanded<-expand_mixture_fit(f)
  ll<-forward_backward_mixture(d$beta_corrected,expanded$initial,expanded$transition,expanded$weights,expanded$means,expanded$sds,lengths)$log_likelihood
  expect_equal(ll,f$log_likelihood,tolerance=1e-9)
  expect_gte(g$log_likelihood+1e-6,f$log_likelihood)
 }
 expect_equal(readfit('G5')$n_parameters-readfit('M221')$n_parameters,14L)
})
test_that('prediction is causal and sequence recursion matches direct enumeration',{
 y<-c(.1,.4,.9,.2,.7,.5);st<-c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE)
 initial<-c(.6,.4);tr<-matrix(c(.8,.2,.3,.7),2,byrow=TRUE)
 emission<-cbind(dnorm(y,.25,.2),dnorm(y,.7,.25))
 a<-fb_sequences_cpp(emission,initial,tr,st)
 direct<-function(ids) {
  paths<-as.matrix(expand.grid(rep(list(1:2),length(ids))))
  sum(apply(paths,1,function(s)initial[s[1]]*prod(tr[cbind(head(s,-1),tail(s,-1))])*prod(emission[cbind(ids,s)])))
 }
 expect_equal(a$log_likelihood,log(direct(1:3))+log(direct(4:6)),tolerance=1e-12)
 changed<-emission;changed[3,]<-rev(changed[3,])
 b<-fb_sequences_cpp(changed,initial,tr,st)
 expect_equal(a$filtered[1:2,],b$filtered[1:2,],tolerance=1e-12)
 expect_equal(a$filtered[4:6,],b$filtered[4:6,],tolerance=1e-12)
 expect_equal(rowSums(a$posterior),rep(1,6),tolerance=1e-12)
 expect_equal(sum(a$transition_sum),4,tolerance=1e-12)
})
test_that('published sequencing extraction retains coordinate and ratio values',{
 ref<-fread(file.path(p,'data/genomics/Chiang2009_sequencing_segments_hg18.csv'))
 expect_equal(nrow(ref[cell_line=='HCC1143']),407L)
 z<-ref[cell_line=='HCC1143' & chromosome==1 & start==210542048]
 expect_equal(z$end,211135191);expect_equal(z$copy_ratio,2.091)
 z<-ref[cell_line=='HCC1143' & chromosome==1 & start==221670740]
 expect_equal(z$end,233147476);expect_equal(z$copy_ratio,1.072)
})
cat('Genomic result checks completed.\n')
