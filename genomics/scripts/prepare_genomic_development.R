#!/usr/bin/env Rscript
root<-normalizePath('genomics')
.libPaths(c(file.path(root,'library'),.libPaths()))
library(data.table)
library(aroma.light)
library(ggplot2)
p<-file.path(root,'data/genomics')
coords<-as.data.table(readRDS(file.path(p,'snp6_coordinates_hg18.rds')))
tumor<-as.data.table(readRDS(file.path(p,'GSM337641_signals.rds')))
normal<-as.data.table(readRDS(file.path(p,'GSM337662_signals.rds')))
stopifnot(nrow(coords)==nrow(tumor),nrow(tumor)==nrow(normal),ncol(tumor)==2,ncol(normal)==2)
print(names(tumor));print(names(normal))
# Aroma's exported PSCN binary stores total signal then allele-B fraction.
stopifnot(grepl('total',names(tumor)[1],ignore.case=TRUE),
          grepl('fracB',names(tumor)[2],ignore.case=TRUE))
d<-copy(coords);d[,marker_index:=.I]
d[,`:=`(tumor_total=tumor[[1]],beta_tumor=tumor[[2]],normal_total=normal[[1]],beta_normal=normal[[2]])]
# Only chromosome 1 is opened for this first development analysis.
d<-d[chromosome==1 & !is.na(chromosome)]
audit<-data.table(stage='chromosome_1_all_units',n=nrow(d))
d<-d[is.finite(position) & position>0 & is.finite(beta_tumor) & is.finite(beta_normal)]
audit<-rbind(audit,data.table(stage='finite_mapped_SNP_measurements',n=nrow(d)))
d<-d[beta_normal>=.3 & beta_normal<=.7]
audit<-rbind(audit,data.table(stage='normal_heterozygote_filter',n=nrow(d)))
setorder(d,position,marker_index)
d<-d[!duplicated(position)]
audit<-rbind(audit,data.table(stage='unique_positions',n=nrow(d)))
d[,beta_corrected:=as.numeric(normalizeTumorBoost(beta_tumor,beta_normal,
                                                 muN=rep(.5,.N),preserveScale=FALSE,flavor='v4'))]
stopifnot(all(is.finite(d$beta_corrected)))
d[,sequence:=cumsum(c(TRUE,diff(position)>1e6))]
d[,total_log_ratio:=log2(tumor_total/normal_total)]
fwrite(d,file.path(p,'development_chr1.csv'))
fwrite(audit,file.path(root,'audit/genomic_development_exclusions.csv'))
fwrite(d[,.(n=.N,start=min(position),end=max(position)),by=sequence],
       file.path(root,'audit/genomic_development_sequences.csv'))
z<-melt(d,id.vars=c('position','marker_index'),measure.vars=c('beta_normal','beta_tumor','beta_corrected','total_log_ratio'))
g<-ggplot(z,aes(position/1e6,value))+geom_point(size=.3,alpha=.45)+
  facet_wrap(~variable,ncol=1,scales='free_y')+theme_bw(base_size=11)+
  labs(x='Chromosome 1 position (Mb, hg18)',y='Measured or corrected signal')
ggsave(file.path(root,'results/genomic_development_raw.png'),g,width=10,height=8,dpi=160)
print(audit);print(summary(d$beta_corrected));print(d[,.(n=.N),by=sequence])
