# six spatial simulators scDesign3 [7], SRTsim [8], scMultiSim [11], stLearn [3], Spider [9], and spaSim [4].
## 1. scDesign3
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("SONGDONGYUAN1994/scDesign3")
library(scDesign3)


file_path = "/Users/oukanyou/Desktop/SpatialSimBench/data_rest/BREAST_sc.rds"
breast_sc_sce <- readRDS(file_path)
breast_sc_sce

example_simu <- scdesign3(
  sce = breast_sc_sce,
  assay_use = "counts",
  celltype = "cellType",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = NULL,
  mu_formula = "s(pseudotime, k = 10, bs = 'cr')",
  sigma_formula = "s(pseudotime, k = 5, bs = 'cr')",
  family_use = "nb",
  n_cores = 2,
  correlation_function = "default",
  usebam = FALSE,
  corr_formula = "1",
  copula = "gaussian",
  fastmvn = FALSE,
  DT = TRUE,
  pseudo_obs = FALSE,
  family_set = c("gauss", "indep"),
  important_feature = "all",
  nonnegative = TRUE,
  return_model = FALSE,
  nonzerovar = FALSE,
  parallelization = "mcmapply",
  BPPARAM = NULL,
  trace = FALSE
)

## 2. SRTsim
devtools::install_github('xzhoulab/SRTsim')
devtools::install_github('xzhoulab/SRTsim')
library(SRTsim)



