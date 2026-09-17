# Compiled equivalent of the public R EM updates; no model or objective change.
library(Rcpp)
sourceCpp('numerics/R/fast_mixture_em.cpp')
fast_em <- function(y, f, counts=f$component_counts, lengths=NULL,
                    max_iter=3000L, tolerance=1e-9, sigma_min=.05) {
  starts <- rep(FALSE,length(y)); starts[cumsum(c(1L,head(if(is.null(lengths)) length(y) else lengths,-1L)))] <- TRUE
  z <- fast_mixture_em_cpp(y,f$initial,f$transition,as.matrix(f$weights),
       as.matrix(f$means),as.matrix(f$sds),as.integer(counts),starts,
       as.integer(max_iter),tolerance,sigma_min)
  z$n_parameters <- mixture_n_parameters(z)
  z
}
