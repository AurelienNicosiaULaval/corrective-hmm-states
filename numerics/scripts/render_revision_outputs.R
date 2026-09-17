#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
source("numerics/scripts/run_elk_application.R")
source("numerics/R/revision_validation.R")
out<-"empirical/revision2"
steps<-read.csv(file.path(out,"steps.csv"))
summary<-read.csv(file.path(out,"model_summary.csv"))
pred<-read.csv(file.path(out,"heldout_predictions.csv"))
parameters<-read.csv(file.path(out,"parameters.csv"))
primary<-summary |> filter(floor==.05,fraction==1)
scores<-pred |> group_by(ID,model) |> summarise(n=n(),log_score=sum(log_score),
  coverage=mean(pit>=.05&pit<=.95),
  KS=max(pmax(seq_along(pit)/n()-sort(pit),sort(pit)-(seq_along(pit)-1)/n())),
  negative_mass=mean(negative_mass),.groups="drop")
folds<-pred |> group_by(ID,fraction,model) |> summarise(n=n(),log_score=sum(log_score),.groups="drop")
write.csv(scores,file.path(out,"heldout_summary.csv"),row.names=FALSE)
write.csv(folds,file.path(out,"heldout_by_fold.csv"),row.names=FALSE)
get_fit<-function(ID,model,floor=.05) {
 n<-sum(steps$ID==ID)
 readRDS(sprintf("%s/fits/elk%s_sd%s_train%s_%s.rds",out,ID,floor,n,model))
}
save_figure<-function(plot,name,width,height) {
 ggsave(file.path(out,paste0(name,".pdf")),plot,width=width,height=height,device=grDevices::cairo_pdf)
}
theme_set(theme_classic(base_size=11)+theme(strip.background=element_blank(),
  strip.text=element_text(hjust=0),legend.position="bottom",plot.title=element_text(size=11)))
state_colors<-c("1"="#286E9F","2"="#BF632F","3"="#72518C","4"="#398774")
model_colors<-c(G3="#555555",M21="#276D9B",M22="#BD622B")

# Every animal is shown, including unfavorable held-out comparisons.
density_rows<-lapply(unique(steps$ID),function(ID) {
 f<-get_fit(ID,"M22"); grid<-seq(-.15,max(steps$log1p_step_km[steps$ID==ID])+.35,length.out=1500)
 cbind(ID,density_decomposition(f,grid))
}) |> bind_rows() |> mutate(state=factor(state),component=factor(component))
total_density<-density_rows |> group_by(ID,y) |> summarise(density=sum(density),.groups="drop")
p_density<-ggplot(steps,aes(log1p_step_km))+
 geom_histogram(aes(y=after_stat(density)),binwidth=.18,boundary=0,fill="grey87",color="white",linewidth=.2)+
 geom_line(data=density_rows,aes(x=y,y=density,color=state,linetype=component),linewidth=.5)+
 geom_line(data=total_density,aes(x=y,y=density),linewidth=.7)+
 facet_wrap(~ID,ncol=1,scales="free_y")+scale_color_manual(values=state_colors)+
 labs(x="log(1 + step length in km)",y="Density",title="A  Four components in two aggregate states",color="State",linetype="Component")
p_pit<-ggplot(pred |> filter(model %in% names(model_colors)),aes(pit,color=model))+
 stat_ecdf(geom="step",linewidth=.6)+geom_abline(slope=1,intercept=0,linetype=2,color="grey55")+
 facet_wrap(~ID,ncol=1)+scale_color_manual(values=model_colors)+coord_equal()+
 labs(x="Held-out one-step PIT",y="Empirical distribution",title="B  Fixed-parameter sequential predictions",color="Model")
save_figure(p_density+p_pit,"elk_observation_diagnostics",8.5,10.0)

# Worked case: three atoms represented by two versus three Markov states.
ID<-"elk-163"; tr<-steps |> filter(.data$ID==.env$ID); f<-get_fit(ID,"M21"); g<-get_fit(ID,"G3")
parts<-density_decomposition(f,seq(-.15,3.6,length.out=1500)) |> mutate(state=factor(state),component=factor(component))
stated<-parts |> group_by(y,state) |> summarise(density=sum(density),.groups="drop")
p_case_density<-ggplot(tr,aes(log1p_step_km))+
 geom_histogram(aes(y=after_stat(density)),binwidth=.12,boundary=0,fill="grey87",color="white")+
 geom_line(data=parts,aes(x=y,y=density,color=state,group=interaction(state,component)),linetype=2)+
 geom_line(data=stated,aes(x=y,y=density,color=state),linewidth=.7)+
 scale_color_manual(values=state_colors)+labs(x="log(1 + step length in km)",y="Density",color="Aggregate state",
 title="A  Elk 163: two components within the lower-movement state")
paths<-bind_rows(data.frame(index=tr$step_index,step=tr$step_km,state=factor(f$viterbi),model="M21: two aggregate states"),
 data.frame(index=tr$step_index,step=tr$step_km,state=factor(g$viterbi),model="G3: three Gaussian states"))
p_paths<-ggplot(paths,aes(index,step))+geom_line(color="grey70",linewidth=.35)+
 geom_point(aes(color=state),size=1.15)+facet_wrap(~model,ncol=1)+scale_color_manual(values=state_colors)+
 labs(x="Consecutive observation interval",y="Step length (km)",color="State",title="B  Corresponding movement sequences")
mapping<-as.data.frame(table(Gaussian=g$viterbi,Aggregate=f$viterbi)) |> group_by(Gaussian) |> mutate(proportion=Freq/sum(Freq))
p_map<-ggplot(mapping,aes(Aggregate,Gaussian,fill=proportion))+geom_tile(color="white")+
 geom_text(aes(label=sprintf("%.2f",proportion)),size=3.3)+scale_fill_gradient(low="white",high="#7DB2D2",limits=c(0,1))+
 labs(x="M21 aggregate state",y="G3 state",fill="Row fraction",title="C  Allocation of the three Gaussian states")
save_figure((p_case_density/p_paths/p_map)+plot_layout(heights=c(1,1.25,.65)),"elk_state_mapping",7.5,9.8)

# Preserve the original all-animal state-mapping comparison in the supplement.
allmapping<-lapply(unique(steps$ID),function(ID) {
 f<-get_fit(ID,"M22"); g<-get_fit(ID,"G3")
 as.data.frame(table(Gaussian=g$viterbi,Aggregate=f$viterbi)) |> mutate(ID=ID) |>
 group_by(ID,Gaussian) |> mutate(proportion=Freq/sum(Freq))
}) |> bind_rows()
save_figure(ggplot(allmapping,aes(Aggregate,Gaussian,fill=proportion))+geom_tile(color="white")+
 geom_text(aes(label=sprintf("%.2f",proportion)),size=3.4)+facet_wrap(~ID,ncol=2)+
 scale_fill_gradient(low="white",high="#7DB2D2",limits=c(0,1))+
 labs(x="M22 aggregate state",y="G3 state",fill="Row fraction"),"elk_all_mappings",7,4.6)

# The four-peak criticism is also shown directly against the variance bound.
sensitivity<-lapply(c(.025,.05,.1),function(floor) {
 fit<-get_fit("elk-115","M22",floor)
 density_decomposition(fit,seq(-.15,3.2,length.out=1500)) |> mutate(floor=floor)
}) |> bind_rows()
sensitivity_total<-sensitivity |> group_by(floor,y) |> summarise(density=sum(density),.groups="drop")
save_figure(ggplot(sensitivity,aes(y,density,color=factor(state),group=interaction(state,component)))+
 geom_line(linetype=2,linewidth=.5)+geom_line(data=sensitivity_total,aes(y,density,group=1),inherit.aes=FALSE)+
 facet_wrap(~floor,ncol=1,labeller=label_both)+scale_color_manual(values=state_colors)+
 labs(x="log(1 + step length in km)",y="Density",color="State"),"elk115_modes_sensitivity",6.8,7.0)

# Produce tables from the same machine-readable results used by the figures.
table_tex<-function(data,path,caption,label,align=NULL,long=FALSE) {
 if(is.null(align))align<-paste0("l",paste(rep("r",ncol(data)-1),collapse=""))
 lines<-apply(data,1,function(x)paste0(paste(x,collapse=" & ")," \\\\"))
 if(long)text<-c(paste0("\\begin{longtable}{",align,"}"),paste0("\\caption{",caption,"}\\label{",label,"}\\\\"),
  "\\toprule",paste0(paste(names(data),collapse=" & ")," \\\\\\midrule\\endfirsthead"),
  "\\toprule",paste0(paste(names(data),collapse=" & ")," \\\\\\midrule\\endhead"),lines,"\\bottomrule\\end{longtable}") else
 text<-c("\\begin{table}[tbp]","\\centering\\small",paste0("\\caption{",caption,"}\\label{",label,"}"),
  paste0("\\begin{tabular}{",align,"}\\toprule"),paste0(paste(names(data),collapse=" & ")," \\\\\\midrule"),
  lines,"\\bottomrule\\end{tabular}","\\end{table}")
 writeLines(text,path)
}
dir.create("article/tables",recursive=TRUE,showWarnings=FALSE)
wide<-primary |> select(ID,model,BIC) |> pivot_wider(names_from=model,values_from=BIC) |> select(ID,G2,G3,G4,M21,M22)
wide$ID<-sub("elk-","",wide$ID)
wide[-1]<-lapply(wide[-1],function(x)sprintf("%.1f",x));names(wide)[1]<-"Elk"
table_tex(wide,"article/tables/elk_model_comparison.tex",
 "BIC for the five models at the primary standard-deviation bound 0.05. G2, G3 and G4 have 7, 14 and 23 parameters; M21 and M22 have 10 and 13. Lower values are preferred.","tab:elk-comparison")
predictive_table<-scores |> mutate(ID=sub("elk-","",ID),log_score=sprintf("%.2f",log_score),
 coverage=sprintf("%.3f",coverage),KS=sprintf("%.3f",KS),negative_mass=sprintf("%.3f",negative_mass))
names(predictive_table)<-c("Elk","Model","$n_{\\rm test}$","Log score","Coverage","PIT distance","$P(Y<0)$")
table_tex(predictive_table,"supplement/tables/elk_predictive_diagnostics.tex",
 "Held-out predictive summaries over two disjoint blocks: observations after the first 60\\% through 80\\%, and after 80\\% through the end. Each block uses parameters fitted to its preceding prefix. Coverage is the fraction of PIT values in $[0.05,0.95]$; PIT distance is the empirical Kolmogorov distance. $P(Y<0)$ is the mean predictive probability assigned below the support of the transformed response. Higher log scores are better.","tab:si-predictions",align="llrrrrr")
param_table<-parameters |> filter(fraction==1,floor==.05,model %in% c("M21","M22")) |>
 transmute(Elk=sub("elk-","",ID),Model=model,State=state,Component=component,
   Mean=sprintf("%.3f",mean),SD=sprintf("%.3f",sd),Weight=sprintf("%.3f",weight),
   Occupation=sprintf("%.3f",occupation),Self=sprintf("%.3f",self_transition))
table_tex(param_table,"supplement/tables/elk_mixture_parameters.tex",
 "Component parameters at bound 0.05. Means and standard deviations are on the transformed scale. Occupation is the average smoothed state probability and Self is the state self-transition probability. Components within a state share the same transition dynamics.","tab:si-parameters",align="llrrrrrrr",long=TRUE)
sen_table<-summary |> filter(fraction==1) |> select(ID,floor,model,BIC) |> pivot_wider(names_from=model,values_from=BIC) |>
 transmute(Elk=sub("elk-","",ID),Bound=sprintf("%.3f",floor),
  d21=sprintf("%.1f",M21-G3),d22=sprintf("%.1f",M22-G3),d24=sprintf("%.1f",M22-G4))
names(sen_table)[3:5]<-c("M21 $-$ G3","M22 $-$ G3","M22 $-$ G4")
table_tex(sen_table,"supplement/tables/elk_sensitivity.tex",
 "BIC differences across standard-deviation bounds. Negative entries favor the mixture model. The first and third contrasts hold the total number of Gaussian components fixed.","tab:si-sensitivity")
con_table<-primary |> transmute(Elk=sub("elk-","",ID),Model=model,
  Gain=format(relative_gain,scientific=TRUE,digits=2),AtBound=floor_hits,
  MinTrans=sprintf("%.4f",min_transition),MaxSelf=sprintf("%.4f",max_self),
  NearBest=paste0(refined_within_001,"/",refined_total))
table_tex(con_table,"supplement/tables/elk_convergence.tex",
 "Numerical audit at bound 0.05. Gain is the final relative likelihood increment; AtBound counts component standard deviations at the bound. MinTrans and MaxSelf summarize transition boundaries. NearBest counts original refined starts within 0.01 log-likelihood units of the final selected solution, including subsequent warm-start improvements. Convergence of EM does not establish a global maximum.","tab:si-convergence",align="llrrrrr")
file.copy(file.path(out,c("elk_observation_diagnostics.pdf","elk_state_mapping.pdf")),"article/figures",overwrite=TRUE)
cat("Tables and figures generated from validated result files.\n")
