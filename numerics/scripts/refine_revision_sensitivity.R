#!/usr/bin/env Rscript
# Cross-floor warm starts prevent artificial reversals from missed local maxima.
library(Rcpp)
source("numerics/scripts/run_elk_application.R")
source("numerics/R/revision_validation.R")
Rcpp::sourceCpp("numerics/R/forward_backward.cpp")
forward_backward_core <- forward_backward_compiled
steps <- prepare_steps(load_elk_data())
models <- list(M21=c(2L,1L),M22=c(2L,2L),G2=c(1L,1L),G3=rep(1L,3),G4=rep(1L,4))
floors <- c(.025,.05,.1)
path_for <- function(ID,n,floor,model)
  sprintf("empirical/revision2/fits/elk%s_sd%s_train%s_%s.rds",ID,floor,n,model)
for(pass in 1:3) for(ID in unique(steps$ID)) {
  y <- steps$log1p_step_km[steps$ID==ID]; n<-length(y)
  for(model in names(models)) for(floor in floors) {
    file <- path_for(ID,n,floor,model); fit <- readRDS(file)
    sources <- lapply(setdiff(floors,floor),function(s)readRDS(path_for(ID,n,s,model)))
    inits <- lapply(sources,canonical_parameters,counts=models[[model]])
    if(model %in% c("G3","G4")) {
      mixture <- if(model=="G3")"M21" else "M22"
      inits <- c(inits,list(expand_mixture_fit(readRDS(path_for(ID,n,floor,mixture)))))
    }
    before <- fit$log_likelihood
    for(i in seq_along(inits)) {
      init <- inits[[i]]; init$sds <- pmax(init$sds,floor)
      candidate <- refine_candidate(y,init,models[[model]],floor,10000L,1e-8)
      fit <- adopt_candidate(fit,candidate,models[[model]],y,paste0("warm_pass",pass,"_",i))
    }
    saveRDS(fit,file)
    if(fit$log_likelihood>before+1e-5)cat(ID,model,floor,"gain",fit$log_likelihood-before,"\n")
  }
}
# Check nested likelihood inequalities and nondecrease as the floor relaxes.
for(ID in unique(steps$ID)) {
  n<-sum(steps$ID==ID)
  for(model in names(models)) {
    ll<-vapply(floors,function(f)readRDS(path_for(ID,n,f,model))$log_likelihood,numeric(1))
    stopifnot(all(diff(ll)<1e-4))
  }
  for(f in floors)for(pair in list(c("M21","G3"),c("M22","G4"))) {
    ll<-vapply(pair,function(m)readRDS(path_for(ID,n,f,m))$log_likelihood,numeric(1))
    stopifnot(diff(ll)>-1e-4)
  }
}
cat("Sensitivity nesting checks passed.\n")
