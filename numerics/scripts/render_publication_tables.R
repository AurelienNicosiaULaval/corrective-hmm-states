# Recreate the two editorial tables from saved numerical outputs.
library(data.table)
x<-fread('genomics/results/genomic_sensitivity/summary.csv')
labels<-c(floor_005='$\\sigma_{\\min}=0.005$',floor_020='$\\sigma_{\\min}=0.020$',thin_2='Every second marker',thin_5='Every fifth marker',normal_025='Normal fraction $[0.25,0.75]$',normal_035='Normal fraction $[0.35,0.65]$')
lines<-c('\\begin{table}[htbp]','\\centering','\\small','\\caption{M221 sensitivity. Agreement is evaluated at markers shared with the primary analysis; parameters are refitted within each scenario.}','\\label{tab:si-genomic-sensitivity}','\\begin{tabular}{lrrr}','\\toprule','Scenario & $n$ & Agreement (\\%) & Runs \\\\','\\midrule')
for(i in seq_len(nrow(x)))lines<-c(lines,sprintf('%s & %d & %.2f & %d \\\\',labels[x$setting[i]],x$n[i],100*x$state_agreement[i],x$runs[i]))
writeLines(c(lines,'\\bottomrule','\\end{tabular}','\\end{table}'),'supplement/tables/genomic_sensitivity.tex')
x<-rbindlist(list(fread('genomics/results/genomic_pilot/metrics.csv')[,.(model,parameters,BIC)],fread('genomics/results/genomic_symmetry_benchmark/metrics.csv')[,.(model,parameters,BIC)]))
y<-fread('genomics/results/genomic_validation/scores.csv')
ord<-c('G2','M21','G3','M221','G5','F2','F3');x<-x[match(ord,model)]
lines<-c('\\begin{table}[htbp]','\\centering','\\small','\\caption{Chromosome-1 estimates and fixed-parameter held-out log scores. Lower BIC and higher log scores are preferred within their respective columns. The Replicate column uses the second array pair on chromosome 1, a technical repeat of the same cell line; it is not biological replication. The chromosome blocks concern that cell line. F2/F3 impose reflection symmetry.}','\\label{tab:genomic-comparison}','\\begin{tabular}{lrrrrr}','\\toprule','Model & $p$ & BIC & Replicate & Chr. 12--17 & Chr. 18--22 \\\\','\\midrule')
for(i in seq_len(nrow(x))){m<-x$model[i];scores<-vapply(c('replication','validation','confirmation'),function(d)y[model==m & dataset==d,log_score],numeric(1));lines<-c(lines,sprintf('%s & %d & %.1f & %.2f & %.2f & %.2f \\\\',m,x$parameters[i],x$BIC[i],scores[1],scores[2],scores[3]))}
writeLines(c(lines,'\\bottomrule','\\end{tabular}','\\end{table}'),'article/tables/genomic_comparison.tex')
cat('Publication tables regenerated from saved results.\n')
