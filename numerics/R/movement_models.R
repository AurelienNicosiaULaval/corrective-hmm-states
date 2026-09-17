# Binned-rate HMMs. Source after compiling movement_likelihood.cpp.
movement_families <- c(gamma=0L,weibull=1L,lognormal=2L,truncated_normal=3L)
movement_logit <- function(probability) {
  if(length(probability)==1L)return(numeric())
  log(pmax(probability[-length(probability)],1e-12)/max(tail(probability,1),1e-12))
}
movement_bounds <- function(counts,family,joint=FALSE,shared_angle=FALSE,rho_max=.95) {
  K<-length(counts);R<-sum(counts)
  n_logits<-(K-1)+K*(K-1)+sum(counts-1)
  if(family=="gamma")component<-cbind(c(-9,log(.05)),c(5,log(200)))
  if(family=="weibull")component<-cbind(c(-9,log(.1)),c(5,log(30)))
  if(family=="lognormal")component<-cbind(c(-9,log(.05)),c(5,log(4)))
  if(family=="truncated_normal")component<-cbind(c(-3,log(.05)),c(5,log(4)))
  low<-c(rep(-12,n_logits),rep(component[,1],R))
  high<-c(rep(12,n_logits),rep(component[,2],R))
  if(joint){
    angular_count<-if(shared_angle)K else R
    low<-c(low,rep(c(-2*pi,0),angular_count));high<-c(high,rep(c(2*pi,rho_max),angular_count))
  }
  list(lower=low,upper=high)
}
movement_initial <- function(data,counts,family,joint,shared_angle=FALSE) {
  K<-length(counts);R<-sum(counts)
  initial<-movement_logit(rep(1/K,K))
  tr<-matrix(rgamma(K*K,1.5),K,K)+diag(runif(K,3,18))
  tr<-tr/rowSums(tr)
  transition<-unlist(lapply(seq_len(K),function(k)movement_logit(tr[k,])))
  weight<-unlist(lapply(counts,function(n)movement_logit(rep(1/n,n))))
  positive<-data$rate[is.finite(data$rate)&data$rate>0]
  centres<-as.numeric(quantile(positive,probs=seq(.05,.9,length.out=R),names=FALSE))*exp(rnorm(R,0,.3))
  if(runif(1)<.5)centres<-sample(centres)
  if(family=="gamma")component<-as.vector(rbind(log(centres),log(runif(R,.5,3))))
  if(family=="weibull")component<-as.vector(rbind(log(centres),log(runif(R,.5,2))))
  if(family=="lognormal")component<-as.vector(rbind(log(centres),log(runif(R,.4,1.5))))
  if(family=="truncated_normal")component<-as.vector(rbind(log1p(centres),log(runif(R,.15,.8))))
  angular<-if(joint)as.vector(rbind(sample(c(0,pi),if(shared_angle)K else R,replace=TRUE),runif(if(shared_angle)K else R,.05,.6))) else numeric()
  c(initial,transition,weight,component,angular)
}
movement_evaluate <- function(par,data,counts,family,joint=FALSE,detail=FALSE,shared_angle=FALSE) {
  movement_model_cpp(par,data$lower,data$upper,data$angle,as.integer(counts),
                     movement_families[[family]],joint,detail,shared_angle)
}
movement_expansion <- function(fit,data) {
  d<-movement_evaluate(fit$par,data,fit$counts,fit$family,fit$joint,TRUE,fit$shared_angle)
  group<-d$group;R<-length(group)
  initial<-d$initial[group]*d$weight
  transition<-d$transition[group,group,drop=FALSE]*matrix(d$weight,R,R,byrow=TRUE)
  component<-as.vector(rbind(d$location,d$scale))
  angular<-if(fit$joint)as.vector(rbind(d$direction,d$rho))else numeric()
  c(movement_logit(initial),unlist(lapply(seq_len(R),function(k)movement_logit(transition[k,]))),component,angular)
}
fit_movement <- function(data,counts,family,seed,n_starts=40L,joint=FALSE,
                         shared_angle=FALSE,nested=NULL,rho_max=.95,extra_starts=list()) {
  set.seed(seed);bounds<-movement_bounds(counts,family,joint,shared_angle,rho_max)
  objective<-function(par)movement_evaluate(par,data,counts,family,joint,FALSE,shared_angle)$nll
  starts<-lapply(seq_len(n_starts),function(i)movement_initial(data,counts,family,joint,shared_angle))
  starts<-c(starts,extra_starts)
  if(!is.null(nested))starts<-c(starts,list(movement_expansion(nested,data)))
  optimize_one<-function(par,iter,eval,tolerance=1e-7) {
    par<-pmax(bounds$lower+1e-9,pmin(bounds$upper-1e-9,par))
    tryCatch(nlminb(par,objective,lower=bounds$lower,upper=bounds$upper,
                    control=list(iter.max=iter,eval.max=eval,rel.tol=tolerance,x.tol=1e-8)),
             error=function(e)list(objective=Inf,par=par,convergence=9L,message=conditionMessage(e),iterations=0L,evaluations=c(0,0)))
  }
  screen<-lapply(starts,optimize_one,iter=70L,eval=180L,tolerance=1e-4)
  rank<-order(vapply(screen,`[[`,numeric(1),"objective"))
  chosen<-head(rank,8)
  if(!is.null(nested))chosen<-unique(c(chosen,length(starts)))
  refined<-lapply(screen[chosen],function(s)optimize_one(s$par,1200L,3500L,1e-9))
  best<-refined[[which.min(vapply(refined,`[[`,numeric(1),"objective"))]]
  # A numerical restart checks the optimizer's own stopping result.
  final<-optimize_one(best$par,2000L,6000L,1e-10)
  if(final$objective<best$objective)best<-final
  rows<-function(xs,stage,indices)do.call(rbind,lapply(seq_along(xs),function(i)data.frame(stage=stage,start=indices[i],nll=xs[[i]]$objective,convergence=xs[[i]]$convergence,iterations=xs[[i]]$iterations)))
  fit<-list(par=best$par,nll=best$objective,counts=counts,family=family,joint=joint,
            shared_angle=shared_angle,seed=seed,n_starts=n_starts,rho_max=rho_max,
            convergence=best$convergence,message=best$message,
            n_parameters=length(best$par),n_intervals=nrow(data),
            boundary_indices=which(abs(best$par-bounds$lower)<1e-4|abs(best$par-bounds$upper)<1e-4),
            audit=rbind(rows(screen,"screen",seq_along(screen)),rows(refined,"refine",chosen),rows(list(final),"restart",0L)))
  fit$BIC<-2*fit$nll+fit$n_parameters*log(sum(is.finite(data$rate)| (joint&is.finite(data$angle))))
  fit
}

# Reproduce published rounding intervals and retain missing components.
prepare_movement_data <- function(raw) {
  result<-lapply(split(raw,raw$Individual),function(x){
    n<-nrow(x);dx<-diff(x$Easting);dy<-diff(x$Northing)
    step<-sqrt(dx^2+dy^2)/1000;rate<-x[['Daily_movement_rate(km/day)']][seq_len(n-1)]
    angle<-x[['Turning_angle_(radians)']][seq_len(n-1)]
    invalid<-step==0|c(TRUE,head(step,-1)==0)
    angle[invalid]<-NA_real_
    data.frame(ID=x$Individual[1],interval=seq_len(n-1),rate=rate,step_km=step,
               lower=pmax(0,rate-.005),upper=rate+.005,angle=angle,
               distance_open_forest=x[['dist_openfor (meters)']][seq_len(n-1)],
               raw_habitat_code=x$Habitat[seq_len(n-1)])
  })
  do.call(rbind,result)
}

movement_predictions <- function(fit,data) {
  details<-movement_evaluate(fit$par,data,fit$counts,fit$family,fit$joint,TRUE,fit$shared_angle)
  marginal<-data;marginal$angle<-NA_real_
  rates<-movement_evaluate(fit$par,marginal,fit$counts,fit$family,fit$joint,TRUE,fit$shared_angle)
  left<-marginal;left$upper<-left$lower;left$lower[is.finite(left$lower)]<-0
  right<-marginal;right$lower[is.finite(right$lower)]<-0
  lower_cdf<-movement_evaluate(fit$par,left,fit$counts,fit$family,fit$joint,TRUE,fit$shared_angle)
  upper_cdf<-movement_evaluate(fit$par,right,fit$counts,fit$family,fit$joint,TRUE,fit$shared_angle)
  predicted<-details$initial
  result<-matrix(NA_real_,nrow(data),4,dimnames=list(NULL,c('log_score','rate_log_score','cdf_lower','cdf_upper')))
  for(t in seq_len(nrow(data))){
    atom_weight<-predicted[details$group]*details$weight
    rate_score<-log(sum(atom_weight*exp(rates$log_components[t,])))
    result[t,]<-c(details$log_score[t],rate_score,
      sum(atom_weight*exp(lower_cdf$log_components[t,])),
      sum(atom_weight*exp(upper_cdf$log_components[t,])))
    z<-details$log_emission[t,];filtered<-predicted*exp(z-max(z))
    predicted<-as.numeric((filtered/sum(filtered))%*%details$transition)
  }
  result[!is.finite(data$rate),2:4]<-NA_real_
  as.data.frame(result)
}

movement_posteriors <- function(fit,data) {
  d<-movement_evaluate(fit$par,data,fit$counts,fit$family,fit$joint,TRUE,fit$shared_angle)
  maximum<-apply(d$log_emission,1,max)
  e<-exp(d$log_emission-maximum)
  fb<-forward_backward_compiled(e,d$initial,d$transition)
  stopifnot(abs(fb$log_likelihood+sum(maximum)+d$nll)<1e-8)
  N<-nrow(data);K<-length(fit$counts)
  score<-matrix(-Inf,N,K);back<-matrix(1L,N,K)
  score[1,]<-log(d$initial)+d$log_emission[1,]
  if(N>1)for(t in 2:N)for(k in seq_len(K)){
    candidate<-score[t-1,]+log(d$transition[,k])
    back[t,k]<-which.max(candidate);score[t,k]<-max(candidate)+d$log_emission[t,k]
  }
  path<-integer(N);path[N]<-which.max(score[N,])
  if(N>1)for(t in (N-1):1)path[t]<-back[t+1,path[t+1]]
  responsibility<-matrix(NA_real_,N,sum(fit$counts))
  for(a in seq_len(ncol(responsibility))){
    k<-d$group[a]
    responsibility[,a]<-fb$posterior[,k]*exp(log(d$weight[a])+d$log_components[,a]-d$log_emission[,k])
  }
  list(details=d,posterior=fb$posterior,viterbi=path,responsibility=responsibility)
}
