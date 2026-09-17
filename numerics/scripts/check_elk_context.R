#!/usr/bin/env Rscript
# Supplementary variables describe the same GPS fixes; they are not fitted labels.
library(Rcpp)
library(dplyr)
library(digest)
library(jsonlite)
source("numerics/scripts/run_elk_application.R")
source("numerics/R/revision_validation.R")
Rcpp::sourceCpp("numerics/R/forward_backward.cpp")
forward_backward_core <- forward_backward_compiled
input <- tempfile(fileext=".txt")
utils::download.file("https://ndownloader.figshare.com/files/5594180",input,quiet=TRUE)
stopifnot(digest::digest(file=input,algo="sha256")==
  "36d423b2832aecd929f7d539da8779965b93a48fa6aa1321fdb7b194e68db2a5")
context <- read.delim(input,check.names=FALSE)
context <- context[context$Individual!="",]
raw <- load_elk_data()
stopifnot(nrow(context)==735L,all(context$Individual==as.character(raw$ID)),
  isTRUE(all.equal(as.matrix(context[,2:3]),as.matrix(raw[,2:3]),check.attributes=FALSE)))
steps <- prepare_steps(raw)
models <- list(M21=c(2L,1L),M22=c(2L,2L),G2=c(1L,1L),G3=rep(1L,3),G4=rep(1L,4))
output <- agreement <- list()
for(ID in unique(steps$ID)) {
  track<-context[context$Individual==ID,]; n<-nrow(track)-1L
  angle<-track[['Turning_angle_(radians)']][seq_len(n)]
  for(model in names(models)) {
    fit<-readRDS(sprintf("empirical/revision2/fits/elk%s_sd0.05_train%s_%s.rds",ID,n,model))
    for(k in seq_len(fit$K)) {
      w<-fit$posterior[,k]; good<-is.finite(angle)
      output[[length(output)+1L]]<-data.frame(ID=ID,model=model,state=k,
        effective_n=sum(w[good]),mean_cos_turn=weighted.mean(cos(angle[good]),w[good]),
        mean_abs_turn=weighted.mean(abs(angle[good]),w[good]),
        mean_distance_open_forest_m=weighted.mean(track[['dist_openfor (meters)']][seq_len(n)],w,na.rm=TRUE),
        occupation=mean(w),median_step_km=median(steps$step_km[steps$ID==ID][fit$viterbi==k]),
        mean_posterior_certainty=mean(apply(fit$posterior,1,max)))
    }
    for(floor in c(.025,.1)) {
      other<-readRDS(sprintf("empirical/revision2/fits/elk%s_sd%s_train%s_%s.rds",ID,floor,n,model))
      # For K=2, account for the only possible label permutation.
      if(fit$K==2L)agreement[[length(agreement)+1]]<-data.frame(ID=ID,model=model,floor=floor,
        agreement=max(mean(other$viterbi==fit$viterbi),mean(3L-other$viterbi==fit$viterbi)))
    }
  }
}
write.csv(bind_rows(output),"empirical/revision2/context_summary.csv",row.names=FALSE)
write.csv(bind_rows(agreement),"empirical/revision2/decoding_sensitivity.csv",row.names=FALSE)
# Published movement rates account for the original observation intervals.
# A targeted robustness check uses elk 163, the worked mixture example.
track<-context[context$Individual=="elk-163",]
rate<-track[['Daily_movement_rate(km/day)']][seq_len(nrow(track)-1)]
stopifnot(length(rate)==158L,all(is.finite(rate)),all(rate>=0))
y<-log1p(rate); fits<-list(); summary<-list()
for(m in seq_along(models)) {
 model<-names(models)[m]
 file<-paste0("empirical/revision2/fits/published_rate_elk163_",model,".rds")
 if(file.exists(file))fit<-readRDS(file)else {
   nested<-if(model=="G3")fits$M21 else if(model=="G4")fits$M22 else NULL
   fit<-fit_audited(y,models[[m]],.05,20261910L+m,nested)
   saveRDS(fit,file)
 }
 fits[[model]]<-fit
 summary[[m]]<-data.frame(ID="elk-163",model=model,n=length(y),log_likelihood=fit$log_likelihood,
                         BIC=fit_bic(fit,length(y)),converged=fit$converged)
}
write.csv(bind_rows(summary),"empirical/revision2/published_rate_sensitivity.csv",row.names=FALSE)
write_json(list(source="https://doi.org/10.6084/m9.figshare.3523667.v1",
  downloaded_file="https://ndownloader.figshare.com/files/5594180",
  sha256="36d423b2832aecd929f7d539da8779965b93a48fa6aa1321fdb7b194e68db2a5",
  alignment="All 735 IDs and coordinate pairs exactly match moveHMM::elk_data.",
  blank_rows_removed=370,angle_alignment="Angle at the departure location of the outgoing step; first angle per animal unavailable.",
  interpretation="Auxiliary descriptive variables from the same trajectory, not independent behavioral validation."),
  "empirical/revision2/context_provenance.json",pretty=TRUE,auto_unbox=TRUE)
cat("Context and published-rate checks complete.\n")
