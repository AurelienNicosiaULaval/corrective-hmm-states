#!/usr/bin/env Rscript
# Run explicitly to install into an isolated library, without replacing global packages.
lib <- file.path('genomics','library')
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(normalizePath(lib),.libPaths()))
options(repos=c(CRAN='https://cloud.r-project.org'))
for(p in c('remotes','BiocManager'))if(!requireNamespace(p,quietly=TRUE))install.packages(p,lib=lib)
for(p in c('Rcpp','data.table','ggplot2','patchwork','clue','testthat','jsonlite','digest'))
 if(!requireNamespace(p,quietly=TRUE))install.packages(p,lib=lib)
BiocManager::install(c('aroma.light','affxparser','DNAcopy'),version='3.21',lib=lib,ask=FALSE,update=FALSE)
remotes::install_github('HenrikBengtsson/sfit@cd8582c949b22ac64bbf0262e7558e45367a9efc',lib=lib,upgrade='never')
for(p in c('aroma.affymetrix','PSCBS'))if(!requireNamespace(p,quietly=TRUE))install.packages(p,lib=lib)
capture.output(sessionInfo(),file='genomics/audit/reproduction_session_info.txt')
