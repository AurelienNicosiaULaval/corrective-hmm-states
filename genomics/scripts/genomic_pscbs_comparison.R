#!/usr/bin/env Rscript
root<-normalizePath('genomics')
.libPaths(c(file.path(root,'library'),.libPaths()))
library(PSCBS)
library(data.table)
out<-file.path(root,'results/genomic_pscbs');dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-fread(file.path(root,'data/genomics/development_chr1.csv'))
pieces<-split(d,d$sequence)
fits<-lapply(seq_along(pieces),function(i) {
  z<-pieces[[i]]
  f<-segmentByPairedPSCBS(CT=z$tumor_total/z$normal_total*2,
    betaT=z$beta_tumor,betaN=z$beta_normal,muN=rep(.5,nrow(z)),
    chromosome=1L,x=z$position,seed=20260910+i,verbose=FALSE)
  saveRDS(f,file.path(out,paste0('sequence',i,'.rds')))
  tab<-as.data.table(getSegments(f,simplify=FALSE));tab[,sequence:=i];tab
})
fwrite(rbindlist(fits,fill=TRUE),file.path(out,'segments.csv'));print(rbindlist(fits,fill=TRUE))
