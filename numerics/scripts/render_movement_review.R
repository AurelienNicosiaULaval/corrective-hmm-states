#!/usr/bin/env Rscript
library(Rcpp)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
Rcpp::sourceCpp('numerics/R/movement_likelihood.cpp')
Rcpp::sourceCpp('numerics/R/forward_backward.cpp')
source('numerics/R/movement_models.R')
out<-'empirical/application_review'
dir.create('article/tables',recursive=TRUE,showWarnings=FALSE)
dir.create('article/figures',recursive=TRUE,showWarnings=FALSE)
dir.create('supplement/tables',recursive=TRUE,showWarnings=FALSE)
dir.create('supplement/figures',recursive=TRUE,showWarnings=FALSE)
data<-read.csv(file.path(out,'movement_context.csv'))
joint<-read.csv(file.path(out,'joint_summary.csv'))
rate<-read.csv(file.path(out,'rate_summary.csv'))
blocks<-read.csv(file.path(out,'joint_validation_blocks.csv'))
scores<-blocks |> group_by(ID,family,model) |> summarise(score=sum(log_score),n=sum(n),.groups='drop')
table_tex<-function(data,caption,label,path,columns=NULL){
  if(is.null(columns))columns<-paste0('l',paste(rep('r',ncol(data)-1),collapse=''))
  value<-lapply(data,function(x)if(is.numeric(x))formatC(x,format='f',digits=1)else gsub('_',' ',x,fixed=TRUE))
  rows<-apply(as.data.frame(value),1,function(x)paste(x,collapse=' & '))
  writeLines(c('\\begin{table}[tbp]','\\centering\\small',paste0('\\caption{',caption,'}'),
               paste0('\\label{',label,'}'),paste0('\\begin{tabular}{',columns,'}'),'\\toprule',
               paste0(paste(names(data),collapse=' & '),' \\\\'),'\\midrule',paste0(rows,' \\\\'),
               '\\bottomrule','\\end{tabular}','\\end{table}'),path)
}
primary<-joint |> filter(family=='weibull',model%in%c('H2','H3','M21')) |> select(ID,model,BIC) |>
  pivot_wider(names_from=model,values_from=BIC)
gains<-scores |> filter(family=='weibull') |> select(ID,model,score) |> pivot_wider(names_from=model,values_from=score) |>
  transmute(ID,`M21--H3`=M21-H3,`M21--H2`=M21-H2)
tab<-left_join(primary,gains,by='ID') |> mutate(ID=sub('elk-','',ID)) |> select(ID,H2,H3,M21,`M21--H3`,`M21--H2`)
names(tab)[1]<-'Elk'
table_tex(tab,'Joint Weibull--wrapped-Cauchy comparison. The first three columns are full-series BIC (smaller is better). The last two are summed predictive log-score differences on five disjoint temporal blocks (positive favors M21). M21, H2 and H3 have 16, 11 and 20 parameters.','tab:elk-joint','article/tables/elk_joint_comparison.tex')
for(type in c('joint','rate')){
  tab<-(if(type=='joint')joint else rate) |> select(ID,family,model,BIC) |> pivot_wider(names_from=model,values_from=BIC) |>
    mutate(ID=sub('elk-','',ID),family=recode(family,truncated_normal='Trunc. normal',lognormal='Lognormal',gamma='Gamma',weibull='Weibull')) |>
    select(ID,family,H2,H3,H4,M21,M22)
  names(tab)[1:2]<-c('Elk','Family')
  table_tex(tab,paste0('Full-series BIC for ',if(type=='joint')'joint rate and angle'else'rate-only',' models. Every comparison uses the same rounding intervals and missing-data treatment.'),
            paste0('tab:movement-',type),paste0('supplement/tables/movement_',type,'.tex'),'llrrrrr')
}
tab<-scores |> select(ID,family,model,score) |> pivot_wider(names_from=model,values_from=score) |>
  mutate(ID=sub('elk-','',ID),family=recode(family,truncated_normal='Trunc. normal'),`M21--H3`=M21-H3) |>
  select(ID,family,H2,H3,M21,`M21--H3`)
names(tab)[1:2]<-c('Elk','Family')
table_tex(tab,'Joint predictive log scores over five disjoint blocks. Parameters are estimated on the preceding training prefix only. These retrospective comparisons are exploratory.','tab:movement-validation','supplement/tables/movement_validation.tex','llrrrr')
sens<-read.csv(file.path(out,'sensitivity_summary.csv')) |> filter(ID=='elk-287')
tab<-sens |> select(family,scenario,model,BIC) |> pivot_wider(names_from=model,values_from=BIC) |>
  left_join(sens |> filter(model=='M21') |> select(family,scenario,agreement),by=c('family','scenario')) |>
  mutate(agreement=100*agreement) |> select(family,scenario,M21,H2,H3,agreement)
names(tab)<-c('Family','Check','M21','H2','H3','Agreement (\\%)')
table_tex(tab,'Elk 287 sensitivity: BIC comparisons within each row and M21 path agreement with the corresponding uncoarsened family. Coarse 005 and 010 group rates below 0.05 and 0.10 km/day and omit angles involving those small movements. Scores across coarsening scenarios concern different observations and must not be compared.','tab:movement-sensitivity','supplement/tables/movement_sensitivity.tex','llrrrr')

x<-data[data$ID=='elk-287',]
f<-readRDS(file.path(out,'fits','elk-287_weibull_joint_M21.rds'));a<-movement_posteriors(f,x);d<-a$details
h<-readRDS(file.path(out,'fits','elk-287_weibull_joint_H3.rds'));b<-movement_posteriors(h,x)
mean_h<-exp(b$details$location)*gamma(1+exp(-b$details$scale));ranks<-rank(mean_h)
x$aggregate<-factor(a$viterbi,levels=1:2,labels=c('Mobile phase','Localized phase'))
x$ordinary<-factor(ranks[b$viterbi],levels=1:3,labels=c('H3: short','H3: moderate','H3: long'))
theme_set(theme_classic(base_size=10)+theme(plot.title=element_text(size=11),legend.position='bottom'))
col<-c('Mobile phase'='#246C91','Localized phase'='#BA562F')
p1<-ggplot(x,aes(interval,rate))+geom_line(color='grey70',na.rm=TRUE)+geom_point(aes(color=aggregate),size=1.1,na.rm=TRUE)+
  scale_y_continuous(trans='log1p')+scale_color_manual(values=col,name=NULL)+labs(x='Recorded interval',y='Rate (km/day)',title='A. Two aggregate phases, three joint components')
p2<-ggplot(x,aes(interval,rate))+geom_line(color='grey70',na.rm=TRUE)+geom_point(aes(color=ordinary),size=1.1,na.rm=TRUE)+
  scale_y_continuous(trans='log1p')+scale_color_manual(values=c('#458358','#BE862E','#705BA5'),name=NULL)+labs(x='Recorded interval',y='Rate (km/day)',title='B. Three ordinary states')
grid<-seq(0.001,log1p(20),length.out=1000);angles<-seq(-pi,pi,length.out=1000)
curves<-bind_rows(lapply(1:3,function(k)data.frame(y=grid,angle=angles,component=factor(k),state=d$group[k],
  rate_density=d$weight[k]*dweibull(expm1(grid),shape=exp(d$scale[k]),scale=exp(d$location[k]))*exp(grid),
  angle_density=d$weight[k]*(1-d$rho[k]^2)/(2*pi*(1+d$rho[k]^2-2*d$rho[k]*cos(angles-d$direction[k]))))))
curves$state<-factor(curves$state,levels=1:2,labels=levels(x$aggregate))
total<-curves |> group_by(y,angle,state) |> summarise(rate_density=sum(rate_density),angle_density=sum(angle_density),.groups='drop')
p3<-ggplot(curves,aes(y,rate_density,color=state))+geom_line(aes(group=component),linetype=2,linewidth=.45)+
  geom_line(data=total,linewidth=.8)+scale_color_manual(values=col)+scale_x_continuous(breaks=log1p(c(0,.1,.5,1,3,10,20)),labels=c('0','.1','.5','1','3','10','20'))+
  labs(x='Rate (km/day), log(1 + rate) scale',y='Density on transformed scale',title='C. State densities and component contributions')+guides(color='none')
p4<-ggplot(curves,aes(angle,angle_density,color=state))+geom_line(aes(group=component),linetype=2,linewidth=.45)+
  geom_line(data=total,linewidth=.8)+scale_color_manual(values=col)+scale_x_continuous(breaks=c(-pi,0,pi),labels=c('-pi','0','pi'))+
  labs(x='Turning angle (radians)',y='Density',title='D. Direction within each phase')+guides(color='none')
combined<-wrap_plots(p1,p2,wrap_plots(p3,p4,nrow=1),ncol=1,heights=c(1,1,1.1))
ggsave(file.path(out,'elk_joint_phase.pdf'),combined,width=8.5,height=9)
file.copy(file.path(out,'elk_joint_phase.pdf'),'article/figures/elk_joint_phase.pdf',overwrite=TRUE)

context<-read.csv(file.path(out,'elk287_phase_path.csv'))
q1<-ggplot(context,aes(interval,probability_mobile))+geom_line(color=col[1],linewidth=.7)+
  labs(x='Recorded interval',y='Smoothed mobile-phase probability',title='A. Conditional state uncertainty')
q2<-ggplot(context,aes(interval,distance_water))+geom_line(color='grey40')+geom_hline(yintercept=1,linetype=2)+
  labs(x='Recorded interval',y='Distance to water (km)',title='B. Habitat context, unused by the mixture fit')
ggsave(file.path(out,'elk_phase_context.pdf'),q1/q2,width=8,height=6)
file.copy(file.path(out,'elk_phase_context.pdf'),'supplement/figures/elk_phase_context.pdf',overwrite=TRUE)
pit<-read.csv(file.path(out,'joint_validation_with_pit.csv')) |> filter(family=='weibull')
q3<-ggplot(pit,aes(pit,color=model))+stat_ecdf(na.rm=TRUE)+geom_abline(slope=1,intercept=0,linetype=2,color='grey50')+
  facet_wrap(~ID)+scale_color_manual(values=c(H2='#707070',H3='#705BA5',M21='#246C91'))+
  labs(x='Randomized predictive PIT',y='Empirical distribution',color='Model')
ggsave(file.path(out,'movement_pit.pdf'),q3,width=8,height=6)
file.copy(file.path(out,'movement_pit.pdf'),'supplement/figures/movement_pit.pdf',overwrite=TRUE)
cat('Application review figures and tables generated from saved fits and predictions.\n')
