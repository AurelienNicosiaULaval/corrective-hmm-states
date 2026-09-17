#!/usr/bin/env Rscript
library(Rcpp)
library(parallel)
library(dplyr)
Rcpp::sourceCpp('numerics/R/movement_likelihood.cpp')
source('numerics/R/movement_models.R')
out <- 'empirical/application_review'
data <- read.csv(file.path(out, 'movement_data.csv'))
files <- list.files(file.path(out, 'fits'), pattern='[.]rds$', full.names=TRUE)
files <- c(files, list.files(file.path(out, 'validation_fits'), pattern='[.]rds$', full.names=TRUE))
refine <- function(path) {
  fit <- readRDS(path)
  id <- sub('_.*', '', basename(path))
  x <- head(data[data$ID == id, ], fit$n_intervals)
  objective <- function(p) movement_evaluate(p,x,fit$counts,fit$family,fit$joint,FALSE,fit$shared_angle)$nll
  bounds <- movement_bounds(fit$counts,fit$family,fit$joint,fit$shared_angle,fit$rho_max)
  original <- fit$nll
  if(fit$convergence != 0L) {
    first <- optim(fit$par,objective,method='L-BFGS-B',lower=bounds$lower,upper=bounds$upper,
                   control=list(maxit=4000,factr=1e6,pgtol=1e-5,ndeps=rep(1e-5,length(fit$par))))
    second <- nlminb(first$par,objective,lower=bounds$lower,upper=bounds$upper,
                     control=list(iter.max=4000,eval.max=15000,rel.tol=1e-10,x.tol=1e-9))
    candidates <- list(list(par=fit$par,nll=fit$nll,convergence=fit$convergence,message=fit$message),
                       list(par=first$par,nll=first$value,convergence=first$convergence,message=first$message),
                       list(par=second$par,nll=second$objective,convergence=second$convergence,message=second$message))
    # Prefer a successful stopping result only within numerical equality of the minimum.
    minimum <- min(vapply(candidates, `[[`, numeric(1), 'nll'))
    eligible <- which(vapply(candidates,function(z)z$nll <= minimum+1e-7 && z$convergence==0L,logical(1)))
    best <- candidates[[if(length(eligible))tail(eligible,1)else which.min(vapply(candidates,`[[`,numeric(1),'nll'))]]
    fit$refinement_history <- list(original_nll=original,original_convergence=fit$convergence,
                                  lbfgsb=first,nlminb=second)
    fit$par <- best$par;fit$nll <- best$nll;fit$convergence <- best$convergence;fit$message <- best$message
    fit$BIC <- 2*fit$nll+fit$n_parameters*log(sum(is.finite(x$rate)|(fit$joint&is.finite(x$angle))))
    fit$boundary_indices <- which(abs(fit$par-bounds$lower)<1e-4|abs(fit$par-bounds$upper)<1e-4)
  }
  gradient <- vapply(seq_along(fit$par),function(j){
    h <- 1e-5;lo <- hi <- fit$par
    lo[j] <- max(bounds$lower[j],lo[j]-h);hi[j] <- min(bounds$upper[j],hi[j]+h)
    (objective(hi)-objective(lo))/(hi[j]-lo[j])
  },numeric(1))
  projected <- gradient
  projected[fit$par <= bounds$lower+1e-5 & gradient>0] <- 0
  projected[fit$par >= bounds$upper-1e-5 & gradient<0] <- 0
  fit$projected_gradient_max <- max(abs(projected))
  saveRDS(fit,path)
  data.frame(file=basename(path),validation=grepl('validation_fits',path),ID=id,family=fit$family,
             joint=fit$joint,n=fit$n_intervals,improvement=original-fit$nll,convergence=fit$convergence,
             max_projected_gradient=fit$projected_gradient_max,bounds=length(fit$boundary_indices))
}
results <- mclapply(files,refine,mc.cores=4L,mc.preschedule=FALSE)
stopifnot(!any(vapply(results,inherits,logical(1),'try-error')))
write.csv(bind_rows(results),file.path(out,'optimization_audit.csv'),row.names=FALSE)
print(bind_rows(results) |> filter(convergence!=0 | max_projected_gradient>.01))
