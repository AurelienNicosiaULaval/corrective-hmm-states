#!/usr/bin/env Rscript
library(Rcpp)
library(data.table)
library(ggplot2)
library(patchwork)
root <- 'genomics'
sourceCpp(file.path(root,'scripts/gaussian_mixture_hmm.cpp'))
source(file.path(root,'scripts/symmetric_genomic_hmm.R'))
out <- file.path(root,'results/genomic_diagnostics'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- fread(file.path(root,'data/genomics/development_chr1.csv'))
m <- readRDS(file.path(root,'results/genomic_pilot/M221.rds'))
g <- readRDS(file.path(root,'results/genomic_pilot/G5.rds'))
labels <- c('Strong imbalance','Moderate imbalance','Approximately balanced')
palette <- c('#9A3F00','#2369A1','#21845C')
d[,state:=m$viterbi];d[,profile:=factor(state,levels=1:3,labels=labels)]
# Transition comparisons are descriptive, not a chi-square test of exact rank.
alignment <- fread(file.path(root,'results/genomic_validation/emission_alignment.csv'))
pairs <- split(alignment$ordinary_state,alignment$mixture_state)
row_tv <- rbindlist(lapply(seq_along(pairs),function(i) {
 j<-pairs[[i]];if(length(j)!=2)return(NULL)
 data.table(profile=labels[i],first=j[1],second=j[2],total_variation=sum(abs(g$transition[j[1],]-g$transition[j[2],]))/2)
}))
fwrite(row_tv,file.path(out,'transition_row_distances.csv'))
fwrite(data.table(index=1:5,singular_value=svd(g$transition)$d),file.path(out,'transition_singular_values.csv'))
fwrite(as.data.table(g$transition),file.path(out,'G5_transition.csv'))
# Reference comparisons preserve assay distinctions.
seg <- fread(file.path(root,'results/genomic_pscbs/segments.csv'))
d[,pscbs_dh:=NA_real_]
for(i in seq_len(nrow(seg)))d[sequence==seg$sequence[i] & position>=seg$dhStart[i] & position<=seg$dhEnd[i],pscbs_dh:=seg$dhMean[i]]
ref <- fread(file.path(root,'data/genomics/Chiang2009_sequencing_segments_hg18.csv'))[cell_line=='HCC1143' & chromosome==1]
d[,sequencing_ratio:=NA_real_]
for(i in seq_len(nrow(ref)))d[position>=ref$start[i] & position<=ref$end[i],sequencing_ratio:=ref$copy_ratio[i]]
fwrite(d[,.(n=.N,pscbs_matched=sum(is.finite(pscbs_dh)),median_pscbs_dh=median(pscbs_dh,na.rm=TRUE),sequencing_matched=sum(is.finite(sequencing_ratio)),median_sequencing_ratio=median(sequencing_ratio,na.rm=TRUE)),by=profile],file.path(out,'assay_summary_by_profile.csv'))
# A pre-existing sequencing interval supplies the zoom boundary, not HMM labels.
zoom <- c(210542048,221670720)
fwrite(d[position>=zoom[1] & position<=zoom[2],.(n=.N,fraction=.N/nrow(d[position>=zoom[1] & position<=zoom[2]]),median_sequencing_ratio=median(sequencing_ratio,na.rm=TRUE)),by=profile],file.path(out,'balanced_gain_region.csv'))
# Predictive PIT uses only the previous forward filter, never smoothed states.
cal <- list();pit_data <- list()
for(key in c('validation','confirmation')) {
 file <- if(key=='validation')'validation_chr12_17.csv' else 'confirmation_chr18_22.csv'
 z <- fread(file.path(root,'data/genomics',file));start<-c(TRUE,diff(z$sequence)!=0)
 for(name in c('G2','M21','G3','G4','M221','G5','F2','F3')) {
  sym<-startsWith(name,'F')
  f<-readRDS(file.path(root,'results',if(sym)'genomic_symmetry_benchmark' else 'genomic_pilot',paste0(name,'.rds')))
  k<-if(sym)length(f$mu) else f$K
  cdf<-emission<-matrix(0,nrow(z),k)
  for(s in seq_len(k)) {
   if(sym) {
    cdf[,s]<-(pnorm(z$beta_corrected,.5-f$mu[s],f$sd[s])+pnorm(z$beta_corrected,.5+f$mu[s],f$sd[s]))/2
    emission[,s]<-(dnorm(z$beta_corrected,.5-f$mu[s],f$sd[s])+dnorm(z$beta_corrected,.5+f$mu[s],f$sd[s]))/2
   } else for(j in seq_len(f$component_counts[s])) {
    cdf[,s]<-cdf[,s]+f$weights[s,j]*pnorm(z$beta_corrected,f$means[s,j],f$sds[s,j])
    emission[,s]<-emission[,s]+f$weights[s,j]*dnorm(z$beta_corrected,f$means[s,j],f$sds[s,j])
   }
  }
  fb<-fb_sequences_cpp(pmax(emission,1e-300),f$initial,f$transition,start)
  pred<-rbind(f$initial,fb$filtered[-nrow(z),,drop=FALSE]%*%f$transition)
  pred[start,]<-matrix(f$initial,sum(start),k,byrow=TRUE)
  u<-rowSums(pred*cdf);r<-qnorm(pmin(1-1e-12,pmax(1e-12,u)))
  prev<-which(!start);ks<-max(abs(sort(u)-(seq_along(u)-.5)/length(u)))+.5/length(u)
  cal[[length(cal)+1]]<-data.table(dataset=key,model=name,n=length(u),coverage95=mean(u>=.025 & u<=.975),pit_KS_distance=ks,residual_lag1=cor(r[prev],r[prev-1]),absolute_residual_lag1=cor(abs(r[prev]),abs(r[prev-1])))
  pit_data[[length(pit_data)+1]]<-data.table(dataset=key,model=name,pit=u)
 }
}
fwrite(rbindlist(cal),file.path(out,'predictive_calibration.csv'));print(rbindlist(cal));print(row_tv)
base_theme <- theme_bw(base_size=10)+theme(legend.position='bottom',panel.grid.minor=element_blank(),plot.title=element_text(size=11))
a <- ggplot(d,aes(position/1e6,beta_corrected,color=profile))+geom_point(size=.2,alpha=.5)+scale_color_manual(values=palette)+labs(x=NULL,y='B-allele fraction',title='A  Three aggregate profiles (M221)',color=NULL)+base_theme+theme(legend.position='none')
bd<-copy(d);bd[,ordinary_state:=factor(g$viterbi,levels=1:5)]
b <- ggplot(bd,aes(position/1e6,beta_corrected,color=ordinary_state))+geom_point(size=.2,alpha=.5)+scale_color_manual(values=c('#9A3F00','#2369A1','#8CBCE0','#21845C','#F6AA59'))+labs(x='Chromosome 1 position (Mb, hg18)',y='B-allele fraction',title='B  Five ordinary Gaussian states (G5)',color='State')+base_theme+theme(legend.position='none')
y<-seq(-.12,1.12,length.out=1400)
comp<-rbindlist(lapply(1:3,function(k)rbindlist(lapply(seq_len(m$component_counts[k]),function(j)data.table(y=y,density=m$weights[k,j]*dnorm(y,m$means[k,j],m$sds[k,j]),profile=factor(labels[k],levels=labels),component=paste(k,j))))))
total<-comp[,.(density=sum(density)),by=.(y,profile)]
c <- ggplot(total,aes(y,density,color=profile))+geom_line(linewidth=.65)+geom_line(data=comp,aes(group=component),linetype=2,linewidth=.35)+scale_color_manual(values=palette)+labs(x='Tumor B-allele fraction',y='Conditional density',title='C  State densities and weighted within-state components',color=NULL)+base_theme
p <- a/b/c+plot_layout(heights=c(1,1,1))
ggsave(file.path(out,'genomic_profiles.pdf'),p,width=7.1,height=7.7);ggsave(file.path(out,'genomic_profiles.png'),p,width=7.1,height=7.7,dpi=170)
# The zoom is an illustration on the development chromosome, not held-out truth.
z<-d[position>=205e6 & position<=235e6]
a<-ggplot(z,aes(position/1e6,beta_corrected,color=profile))+geom_point(size=.4,alpha=.65)+scale_color_manual(values=palette,drop=FALSE)+labs(x=NULL,y='B-allele fraction',title='A  Allelic profiles',color=NULL)+base_theme+theme(legend.position='none')
sz<-seg[dhEnd>=205e6 & dhStart<=235e6]
b<-ggplot(sz,aes(x=pmax(dhStart,205e6)/1e6,xend=pmin(dhEnd,235e6)/1e6,y=dhMean,yend=dhMean))+geom_segment(color='#574C6B',linewidth=.7)+labs(x=NULL,y='Allelic imbalance (DH)',title='B  PSCBS segmentation of the same arrays')+base_theme
rz<-ref[end>=205e6 & start<=235e6]
c<-ggplot(rz,aes(x=pmax(start,205e6)/1e6,xend=pmin(end,235e6)/1e6,y=copy_ratio,yend=copy_ratio))+geom_segment(color='#333333',linewidth=.7)+geom_hline(yintercept=1,linetype=3)+labs(x='Chromosome 1 position (Mb, hg18)',y='Relative copy ratio',title='C  Published independent sequencing estimates')+base_theme
for(n in c('a','b','c'))assign(n,get(n)+scale_x_continuous(limits=c(205,235),expand=expansion(mult=0)))
p<-a/b/c
ggsave(file.path(out,'genomic_assays.pdf'),p,width=7.1,height=6.3);ggsave(file.path(out,'genomic_assays.png'),p,width=7.1,height=6.3,dpi=170)
pits<-rbindlist(pit_data)[model%in%c('M221','G5','F3')]
p<-ggplot(pits,aes(pit))+geom_histogram(aes(y=after_stat(density)),breaks=seq(0,1,.05),fill='#4B7C9B',color='white')+geom_hline(yintercept=1,linetype=2)+facet_grid(dataset~model)+labs(x='One-step predictive PIT',y='Density')+base_theme
ggsave(file.path(out,'genomic_pit.pdf'),p,width=7.1,height=4.6);ggsave(file.path(out,'genomic_pit.png'),p,width=7.1,height=4.6,dpi=160)
