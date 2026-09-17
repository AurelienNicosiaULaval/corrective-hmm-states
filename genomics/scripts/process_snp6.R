#!/usr/bin/env Rscript
# Process the first two accession-ordered HCC1143 tumor-normal pairs.
# Pair selection precedes inspecting their chromosome profiles or model fits.
root <- normalizePath('genomics')
.libPaths(c(file.path(root,'library'),.libPaths()))
library(aroma.affymetrix)
library(data.table)
options(timeout=300)
dir.create(file.path(root,'genomic_processing'),recursive=TRUE,showWarnings=FALSE)
setwd(file.path(root,'genomic_processing'))
raw <- file.path('rawData','GSE13372','GenomeWideSNP_6')
dir.create(raw,recursive=TRUE,showWarnings=FALSE)
ids <- c('GSM337641','GSM337662','GSM337642','GSM337663')
for(id in ids) {
  dest <- file.path(raw,paste0(id,'.CEL'))
  if(!file.exists(dest)) R.utils::gunzip(file.path(root,'data/genomics',paste0(id,'.CEL.gz')),
                                      destname=dest,remove=FALSE,overwrite=FALSE)
}
cs <- AffymetrixCelSet$byName('GSE13372',chipType='GenomeWideSNP_6,Full')
stopifnot(length(cs)==4L)
print(cs)
cache <- file.path(root,'data/genomics','ascrma_object.rds')
ds <- if(file.exists(cache))readRDS(cache)else doASCRMAv2(cs,verbose=-10)
saveRDS(ds,file.path(root,'data/genomics','ascrma_object.rds'))
print(ds)
stopifnot(identical(getNames(ds$total),getNames(ds$fracB)))
for(i in seq_len(length(ds$total))) {
  f <- ds$total[[i]]
  print(f)
  total <- readDataFrame(f)[[1]]
  fracB <- readDataFrame(ds$fracB[[i]])[[1]]
  stopifnot(length(total)==length(fracB))
  dat <- data.frame(total=total,fracB=fracB)
  saveRDS(dat,file.path(root,'data/genomics',paste0(getName(f),'_signals.rds')))
  str(dat)
}
ugp <- AromaUgpFile$byChipType('GenomeWideSNP_6,Full',tags='na30,hg18')
saveRDS(readDataFrame(ugp),file.path(root,'data/genomics','snp6_coordinates_hg18.rds'))
capture.output(sessionInfo(),file=file.path(root,'audit','genomic_session_info.txt'))
