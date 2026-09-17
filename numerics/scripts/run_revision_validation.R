#!/usr/bin/env Rscript
# Run from reproducibility_package. All starts and fixed-parameter forecasts saved.
library(Rcpp)
library(parallel)
library(dplyr)
library(jsonlite)
source("numerics/scripts/run_elk_application.R")
source("numerics/R/revision_validation.R")
reference_forward_backward <- forward_backward_core
Rcpp::sourceCpp("numerics/R/forward_backward.cpp")
# Verify acceleration before using it in any statistical fit.
set.seed(914)
for (K in 1:4) for (T in c(1L,2L,20L)) {
  emission <- matrix(runif(T*K,.0001,1),T,K)
  initial <- normalize_probability(runif(K))
  transition <- normalize_rows(matrix(runif(K*K),K,K))
  stopifnot(isTRUE(all.equal(reference_forward_backward(emission,initial,transition),
                            forward_backward_compiled(emission,initial,transition),
                            tolerance=1e-12)))
}
forward_backward_core <- forward_backward_compiled
dir.create("empirical/revision2/fits",recursive=TRUE,showWarnings=FALSE)
steps <- prepare_steps(load_elk_data())
utils::write.csv(steps[,c("ID","step_index","step_km","log1p_step_km")],
                 "empirical/revision2/steps.csv",row.names=FALSE)
models <- list(M21=c(2L,1L),M22=c(2L,2L),G2=c(1L,1L),G3=rep(1L,3),G4=rep(1L,4))
jobs <- expand.grid(ID=unique(steps$ID), floor=c(.025,.05,.1), fraction=1,
                    stringsAsFactors=FALSE)
jobs <- rbind(jobs,expand.grid(ID=unique(steps$ID),floor=.05,fraction=c(.6,.8),
                              stringsAsFactors=FALSE))
run_job <- function(j) {
  job <- jobs[j,]; track <- steps[steps$ID==job$ID,]
  n <- nrow(track); n_train <- if(job$fraction==1) n else floor(n*job$fraction)
  y <- track$log1p_step_km[seq_len(n_train)]
  result <- list()
  for (m in seq_along(models)) {
    model <- names(models)[m]
    key <- sprintf("elk%s_sd%s_train%s_%s",job$ID,job$floor,n_train,model)
    file <- file.path("empirical/revision2/fits",paste0(key,".rds"))
    if (file.exists(file)) fit <- readRDS(file) else {
      nested <- if(model=="G3") result$M21 else if(model=="G4") result$M22 else NULL
      fit <- fit_audited(y,models[[m]],job$floor,20260910L+j*100L+m,nested)
      saveRDS(fit,file)
    }
    result[[model]] <- fit
    cat(key,sprintf(" ll=%.4f converged=%s\n",fit$log_likelihood,fit$converged))
  }
  result
}
cores <- min(4L,parallel::detectCores())
results <- parallel::mclapply(seq_len(nrow(jobs)),run_job,mc.cores=cores,mc.preschedule=FALSE)
stopifnot(!any(vapply(results,inherits,logical(1),"try-error")))
summaries <- audits <- predictions <- parameters <- list()
for(j in seq_len(nrow(jobs))) {
  job <- jobs[j,]; track <- steps[steps$ID==job$ID,]; n <- nrow(track)
  n_train <- if(job$fraction==1) n else floor(n*job$fraction)
  for(model in names(models)) {
    fit <- results[[j]][[model]]
    identity <- data.frame(ID=job$ID,floor=job$floor,fraction=job$fraction,
                           n_train=n_train,model=model)
    modes <- count_modes(fit)
    refined <- subset(fit$audit,stage=="refine")
    summaries[[length(summaries)+1]] <- cbind(identity,data.frame(
      log_likelihood=fit$log_likelihood,n_parameters=fit$n_parameters,
      BIC=fit_bic(fit,n_train),converged=fit$converged,relative_gain=relative_gain(fit),
      iterations=fit$n_iter,n_modes=modes$count,
      modes=paste(round(modes$locations,5),collapse=";"),
      floor_hits=sum(fit$sds<=job$floor*(1+1e-6)),
      min_transition=min(fit$transition),max_self=max(diag(fit$transition)),
      refined_within_001=sum(refined$log_likelihood>=fit$log_likelihood-.01),
      refined_total=nrow(refined)))
    audits[[length(audits)+1]] <- cbind(identity,fit$audit)
    pred <- sequential_predictions(track$log1p_step_km,fit)
    if(job$fraction<1) {
      end <- if(job$fraction==.6) floor(.8*n) else n
      idx <- seq.int(n_train+1,end)
      predictions[[length(predictions)+1]] <- cbind(identity,step_index=idx,pred[idx,])
    }
    for(k in seq_len(fit$K)) for(a in seq_len(fit$component_counts[k]))
      parameters[[length(parameters)+1]] <- cbind(identity,state=k,component=a,
        mean=fit$means[k,a],sd=fit$sds[k,a],weight=fit$weights[k,a],
        occupation=mean(fit$posterior[,k]),self_transition=fit$transition[k,k])
  }
}
write.csv(bind_rows(summaries),"empirical/revision2/model_summary.csv",row.names=FALSE)
write.csv(bind_rows(audits),"empirical/revision2/start_audit.csv",row.names=FALSE)
write.csv(bind_rows(predictions),"empirical/revision2/heldout_predictions.csv",row.names=FALSE)
write.csv(bind_rows(parameters),"empirical/revision2/parameters.csv",row.names=FALSE)
write_json(list(date=as.character(Sys.Date()),n_starts=200,screen_iterations=70,
  refine_top=20,tolerance=1e-8,maximum_refinement_iterations=10000,
  floors=c(.025,.05,.1),training_fractions=c(.6,.8),seeds="20260910 + 100*job + model",
  compiled_reference_test="12 cases passed at tolerance 1e-12"),
  "empirical/revision2/protocol.json",pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file="empirical/revision2/sessionInfo.txt")
cat("Revision analyses complete.\n")
