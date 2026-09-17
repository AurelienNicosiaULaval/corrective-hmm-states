library(Rcpp);library(dplyr);library(ggplot2);library(jsonlite)
source('numerics/R/hmm_em.R');source('numerics/R/mixture_hmm_em.R');source('numerics/R/simulate.R')
source('numerics/R/optimization_audit.R')
out<-'numerics/output/optimization_audit';dir.create(out,recursive=TRUE,showWarnings=FALSE)
oldfile<-file.path(out,'original_replications.csv')
if(!file.exists(oldfile))file.copy('numerics/output/review_replications.csv',oldfile)
old<-read.csv(oldfile)
tasks<-old[order(-old$T,old$rep),]
run_one<-function(i){
 task<-tasks[i,];path<-file.path(out,sprintf('fit_%d_%03d.rds',task$T,task$rep))
 if(file.exists(path))return(path)
 y<-simulate_hmm(task$T,task$seed)$y;fs<-fit_simulation_audited(y,task$seed)
 vals<-vapply(fs[1:3],`[[`,numeric(1),'log_likelihood')
 baseline<-as.numeric(task[c('logLik_gauss_K2','logLik_gauss_K3','logLik_mix_K2_M2')])
 if(any(vals<baseline-1e-5))stop('Likelihood regression: ',task$T,'/',task$rep,': ',paste(vals-baseline,collapse=','))
 m<-fs$M22;o<-order(rowSums(m$weights*m$means));b<-vapply(fs[1:3],fit_bic,numeric(1),T=task$T)
 row<-data.frame(T=task$T,rep=task$rep,seed=task$seed,logLik_gauss_K2=vals[1],logLik_gauss_K3=vals[2],logLik_mix_K2_M2=vals[3],BIC_gauss_K2=b[1],BIC_gauss_K3=b[2],BIC_mix_K2_M2=b[3],delta_BIC_gauss_K3_minus_K2=b[2]-b[1],delta_BIC_mix_K2_minus_gauss_K3=b[3]-b[2],gaussian_BIC_favours_K3=b[2]<b[1],enriched_BIC_favours_K2_mix=b[3]<b[2],mix_state1_mean=mixture_state_means(m)[o[1]],mix_state2_mean=mixture_state_means(m)[o[2]],mix_state1_occ=colMeans(m$posterior)[o[1]],mix_state2_occ=colMeans(m$posterior)[o[2]])
 for(j in seq_along(fs))fs[[j]]$posterior<-NULL
 saveRDS(list(row=row,fits=fs),paste0(path,'.tmp'));file.rename(paste0(path,'.tmp'),path)
 cat('Completed ',task$T,'/',task$rep,'\n',sep='');path
}
workers<-as.integer(Sys.getenv('AUDIT_WORKERS','6'))
paths<-parallel::mclapply(seq_len(nrow(tasks)),run_one,mc.cores=workers,mc.preschedule=FALSE)
if(any(vapply(paths,inherits,logical(1),'try-error')))stop('Some fits failed; checkpoints retained.')
r<-do.call(rbind,lapply(paths,function(p)readRDS(p)$row));r<-r[order(r$T,r$rep),];rownames(r)<-NULL
write.csv(r,'numerics/output/review_replications.csv',row.names=FALSE)
s<-r |> group_by(T) |> summarise(n_replications=n(),mean_delta_BIC_gauss_K3_minus_K2=mean(delta_BIC_gauss_K3_minus_K2),prop_gaussian_BIC_favours_K3=mean(gaussian_BIC_favours_K3),mean_delta_BIC_mix_K2_minus_gauss_K3=mean(delta_BIC_mix_K2_minus_gauss_K3),prop_enriched_BIC_favours_K2_mix=mean(enriched_BIC_favours_K2_mix),.groups='drop')
write.csv(s,'numerics/output/review_summary_by_T.csv',row.names=FALSE)
comparison<-merge(old,r,by=c('T','rep','seed'),suffixes=c('_old','_new'));write.csv(comparison,file.path(out,'comparison.csv'),row.names=FALSE)
source('numerics/R/plotting.R')
plot_bic_boxplot(r$delta_BIC_gauss_K3_minus_K2[r$T==1500],'article/figures');plot_bic_by_sample_size(s,'article/figures')
# Load the rendering functions without running the legacy numerical protocol.
p<-parse('numerics/scripts/run_review_simulation.R');eval(p[-length(p)])
write_replication_tables(s)
write_json(list(seed=20260610,replications=600,random_starts=10,refine_top=3,max_iterations_per_chunk=3000,tolerance=1e-9,nested_starts='G2 split into M21; M22 reduced into M21; M21 exactly expanded into G3 and duplicated into M22',oracle_starts=FALSE),file.path(out,'protocol.json'),pretty=TRUE,auto_unbox=TRUE)
print(s);cat('All 600 corrected replications and derived outputs saved.\n')
