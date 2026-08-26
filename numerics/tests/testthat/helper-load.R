ROOT <- normalizePath(file.path("..", "..", ".."), mustWork = TRUE)

source(file.path(ROOT, "numerics/R/hmm_em.R"))
source(file.path(ROOT, "numerics/R/mixture_hmm_em.R"))
source(file.path(ROOT, "numerics/R/simulate.R"))
source(file.path(ROOT, "numerics/R/diagnostics.R"))
