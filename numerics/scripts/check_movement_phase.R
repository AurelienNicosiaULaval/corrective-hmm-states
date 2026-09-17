#!/usr/bin/env Rscript
library(Rcpp)
library(dplyr)
library(tidyr)
library(digest)
Rcpp::sourceCpp('numerics/R/movement_likelihood.cpp')
Rcpp::sourceCpp('numerics/R/forward_backward.cpp')
source('numerics/R/movement_models.R')
out<-'empirical/application_review';data<-read.csv(file.path(out,'movement_context.csv'))
x<-data[data$ID=='elk-287',]
f<-readRDS(file.path(out,'fits','elk-287_weibull_joint_M21.rds'))
a<-movement_posteriors(f,x)
mobile<-which.max(vapply(1:2,function(k)max(a$details$location[a$details$group==k]),numeric(1)))
phase<-ifelse(a$viterbi==mobile,'mobile','localized')
input<-tempfile(fileext='.txt');download.file('https://ndownloader.figshare.com/files/5594180',input,quiet=TRUE)
stopifnot(digest(file=input,algo='sha256')=='36d423b2832aecd929f7d539da8779965b93a48fa6aa1321fdb7b194e68db2a5')
raw<-read.delim(input,check.names=FALSE);raw<-raw[raw$Individual=='elk-287'&raw$Individual!='',]
coordinates<-as.matrix(raw[,c('Easting','Northing')])/1000
context<-bind_rows(lapply(unique(phase),function(label){
  rows<-which(phase==label)
  points<-coordinates[unique(c(rows,rows+1)),,drop=FALSE]
  data.frame(phase=label,start=min(rows),end=max(rows),intervals=length(rows),
    median_rate=median(x$rate[rows],na.rm=TRUE),mean_rate=mean(x$rate[rows],na.rm=TRUE),
    maximum_rate=max(x$rate[rows],na.rm=TRUE),path_length_km=sum(x$step_km[rows]),
    spatial_diameter_km=max(dist(points)),
    median_distance_water_km=median(x$distance_water[rows]),
    fraction_within_1km_water=mean(x$distance_water[rows]<1),
    observed_mean_cos=mean(cos(x$angle[rows]),na.rm=TRUE))
}))
write.csv(context,file.path(out,'elk287_phase_context.csv'),row.names=FALSE)
write.csv(data.frame(interval=x$interval,phase=phase,probability_mobile=a$posterior[,mobile],
                     distance_water=x$distance_water),file.path(out,'elk287_phase_path.csv'),row.names=FALSE)

# Sensitivity to the artificial finite-logit approximation to an absorbing state.
objective<-function(p)movement_evaluate(p,x,f$counts,f$family,TRUE)$nll
rows<-list()
for(limit in c(8,12,16)){
 bounds<-movement_bounds(f$counts,f$family,TRUE)
 n_logits<-4L;bounds$lower[1:n_logits]<- -limit;bounds$upper[1:n_logits]<-limit
 start<-pmax(bounds$lower,pmin(bounds$upper,f$par))
  z<-nlminb(start,objective,lower=bounds$lower,upper=bounds$upper,
            control=list(iter.max=5000,eval.max=20000,rel.tol=1e-11))
 if(z$convergence!=0){
   restart<-optim(z$par,objective,method='L-BFGS-B',lower=bounds$lower,upper=bounds$upper,
                   control=list(maxit=4000,factr=1e6,ndeps=rep(1e-5,length(z$par))))
   if(restart$value<=z$objective+1e-7 && restart$convergence==0){
     z<-list(par=restart$par,objective=restart$value,convergence=restart$convergence,message=restart$message)
   }
 }
 ff<-f;ff$par<-z$par;ff$nll<-z$objective;b<-movement_posteriors(ff,x)
 rows[[as.character(limit)]]<-data.frame(logit_bound=limit,nll=z$objective,convergence=z$convergence,
             probability_return=b$details$transition[3-mobile,mobile],
             path_agreement=mean(b$viterbi==a$viterbi))
 saveRDS(list(fit=ff,optimization=z),file.path(out,paste0('elk287_transition_bound_',limit,'.rds')))
}
write.csv(bind_rows(rows),file.path(out,'elk287_transition_boundary.csv'),row.names=FALSE)

# A common randomization for rounding-aware PITs across every fitted model.
pred<-read.csv(file.path(out,'joint_validation_predictions.csv'))
keys<-distinct(pred,ID,interval);set.seed(202609106L);keys$u<-runif(nrow(keys))
pred<-left_join(pred,keys,by=c('ID','interval')) |> mutate(pit=cdf_lower+u*(cdf_upper-cdf_lower))
summary<-pred |> group_by(ID,family,model) |> summarise(n=sum(is.finite(pit)),
  central_90_fraction=mean(pit>=.05&pit<=.95,na.rm=TRUE),mean_pit=mean(pit,na.rm=TRUE),
  rate_log_score=sum(rate_log_score,na.rm=TRUE),joint_log_score=sum(log_score),.groups='drop')
write.csv(summary,file.path(out,'joint_predictive_calibration.csv'),row.names=FALSE)
write.csv(pred,file.path(out,'joint_validation_with_pit.csv'),row.names=FALSE)
pred$normal_pit <- qnorm(pmin(pmax(pred$pit,1e-7),1-1e-7))
dependence <- pred |> group_by(ID,family,model,fraction) |> arrange(interval,.by_group=TRUE) |>
  mutate(previous_pit=lag(normal_pit)) |> ungroup() |> group_by(ID,family,model) |>
  summarise(pit_residual_lag1=cor(normal_pit,previous_pit,use='complete.obs'),.groups='drop')
write.csv(dependence,file.path(out,'predictive_dependence.csv'),row.names=FALSE)
print(context);print(bind_rows(rows));print(summary |> filter(ID=='elk-287',family=='weibull'))
