setwd("../")

library(sharp)
library(glmnet)

# Loading the arguments
args <- commandArgs(trailingOnly = TRUE)
m <- as.character(args[1])
id <- as.numeric(args[2])
ttd <- as.numeric(args[3])

# Loading the data
dnam <- readRDS("Data/imputed_denoised_450K")
covars <- readRDS("Data/covariates.rds")
print(all(rownames(dnam) == rownames(covars)))

# Removing short ttd
ids <- which((covars$case_lc == 0) | ((covars$case_lc == 1) & (covars$ttd >= ttd)))
covars <- covars[ids, ]
dnam <- dnam[rownames(covars), ]
print(dim(dnam))

# Model preparation
covars$gender <- ifelse(as.character(covars$gender) == "Male", yes = 1, no = 0)
xdata <- cbind(dnam, age = covars$age.recr, sex = covars$gender)
penalty <- c(rep(1, ncol(dnam)), 0, 0)
if (m == 2) {
  xdata <- cbind(xdata, smoking = covars$csi)
  penalty <- c(penalty, 0)
}
print(m)
print(length(penalty))

# Variable selection
print(system.time({
  stab <- VariableSelection(
    xdata = xdata, ydata = covars$case_lc,
    family = "binomial", penalty.factor = penalty,
    Lambda = LambdaSequence(lmax = 0.18, lmin = 1e-15, cardinal = 100),
    seed = id, K = 20, tau = 0.8
  )
}))

# Saving outputs
saveRDS(stab, paste0("Results/Iterations/stability_80_m", m, "_", ttd, "_", id, ".rds"))
