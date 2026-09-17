#!/usr/bin/env Rscript
library(Rcpp)
library(dplyr)
library(ggplot2)
library(patchwork)
Rcpp::sourceCpp('numerics/R/movement_likelihood.cpp')
Rcpp::sourceCpp('numerics/R/forward_backward.cpp')
source('numerics/R/movement_models.R')
out<-'empirical/application_review';data<-read.csv(file.path(out,'movement_data.csv'))
params<-states<-paths<-list()
for(family in names(movement_families))for(joint in c(FALSE,TRUE))for(ID in unique(data$ID))for(model in c('M21','H2','H3')){
 name<-paste(ID,family,if(joint)'joint'else'rate',model,sep='_')
 f<-readRDS(file.path(out,'fits',paste0(name,'.rds')));x<-data[data$ID==ID,]
 a<-movement_posteriors(f,x);d<-a$details
 if(family=='weibull')mean_rate<-exp(d$location)*gamma(1+exp(-d$scale))
 if(family=='gamma')mean_rate<-exp(d$location)
 if(family=='lognormal')mean_rate<-exp(d$location+exp(2*d$scale)/2)
 if(family=='truncated_normal')mean_rate<-exp(d$location+exp(2*d$scale)/2+
    pnorm(0,d$location+exp(2*d$scale),exp(d$scale),lower.tail=FALSE,log.p=TRUE)-
    pnorm(0,d$location,exp(d$scale),lower.tail=FALSE,log.p=TRUE))-1
 params[[length(params)+1]]<-data.frame(ID=ID,family=family,joint=joint,model=model,
   component=seq_along(d$weight),state=d$group,weight=d$weight,location=d$location,log_scale=d$scale,
   mean_rate=mean_rate,direction=atan2(sin(d$direction),cos(d$direction)),rho=d$rho)
 for(k in seq_along(f$counts))states[[length(states)+1]]<-data.frame(ID=ID,family=family,joint=joint,model=model,state=k,
   occupation=mean(a$posterior[,k]),self_transition=d$transition[k,k],
   observed_run_count=sum(diff(c(FALSE,a$viterbi==k))==1),
   weighted_cos_turn=weighted.mean(cos(x$angle),a$posterior[,k],na.rm=TRUE),
   weighted_distance_open_forest=weighted.mean(x$distance_open_forest,a$posterior[,k],na.rm=TRUE))
 paths[[length(paths)+1]]<-data.frame(ID=ID,family=family,joint=joint,model=model,interval=x$interval,
   state=a$viterbi,certainty=apply(a$posterior,1,max))
}
write.csv(bind_rows(params),file.path(out,'component_parameters.csv'),row.names=FALSE)
write.csv(bind_rows(states),file.path(out,'state_descriptions.csv'),row.names=FALSE)
write.csv(bind_rows(paths),file.path(out,'decoded_paths.csv'),row.names=FALSE)

# Diagnostic display, kept separate from the manuscript pending assessment.
ID<-'elk-287';x<-data[data$ID==ID,]
f<-readRDS(file.path(out,'fits','elk-287_weibull_joint_M21.rds'));a<-movement_posteriors(f,x)
z<-bind_rows(lapply(1:2,function(k)data.frame(interval=x$interval,state=factor(k),probability=a$posterior[,k])))
theme_set(theme_classic(base_size=11))
colors<-c('1'='#286E9F','2'='#BF632F')
p1<-ggplot(x,aes(interval,rate))+geom_line(color='gray65')+geom_point(aes(color=factor(a$viterbi)),size=1.2)+
 scale_y_continuous(trans='log1p')+scale_color_manual(values=colors,name='Aggregate state')+labs(x='Observation interval',y='Published rate (km/day)',title='Elk 287: rates and aggregate states')
p2<-ggplot(z,aes(interval,probability,color=state))+geom_line()+scale_color_manual(values=colors)+
 labs(x='Observation interval',y='Smoothed probability',title='One dominant phase transition')+guides(color='none')
ggsave(file.path(out,'elk287_phase_diagnostic.pdf'),p1/p2,width=9,height=7)
cat('Movement descriptions saved.\n')
