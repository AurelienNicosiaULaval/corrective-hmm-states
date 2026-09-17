#!/usr/bin/env Rscript
root<-normalizePath('genomics')
.libPaths(c(file.path(root,'library'),.libPaths()))
library(data.table)
library(aroma.light)
p<-file.path(root,'data/genomics')
coords<-as.data.table(readRDS(file.path(p,'snp6_coordinates_hg18.rds')))
prepare<-function(tid,nid,chromosomes,name) {
  tumor<-readRDS(file.path(p,paste0(tid,'_signals.rds')))
  normal<-readRDS(file.path(p,paste0(nid,'_signals.rds')))
  stopifnot(nrow(tumor)==nrow(coords),nrow(normal)==nrow(coords))
  d<-copy(coords);d[,marker_index:=.I]
  d[,`:=`(tumor_total=tumor$total,beta_tumor=tumor$fracB,normal_total=normal$total,beta_normal=normal$fracB)]
  d<-d[chromosome%in%chromosomes & is.finite(position) & position>0 &
         is.finite(beta_tumor) & is.finite(beta_normal) & beta_normal>=.3 & beta_normal<=.7]
  setorder(d,chromosome,position,marker_index);d<-d[!duplicated(d[,.(chromosome,position)])]
  d[,beta_corrected:=as.numeric(normalizeTumorBoost(beta_tumor,beta_normal,
    muN=rep(.5,.N),preserveScale=FALSE,flavor='v4'))]
  stopifnot(all(is.finite(d$beta_corrected)))
  d[,sequence:=cumsum(c(TRUE,diff(chromosome)!=0 | diff(position)>1e6))]
  d[,total_log_ratio:=log2(tumor_total/normal_total)]
  fwrite(d,file.path(p,paste0(name,'.csv')))
  print(d[,.(n=.N),by=chromosome])
}
prepare('GSM337641','GSM337662',18:22,'confirmation_chr18_22')
