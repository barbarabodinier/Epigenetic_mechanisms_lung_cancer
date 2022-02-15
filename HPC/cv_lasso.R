setwd("../")

library(glmnet)

# Loading the arguments
args <- commandArgs(trailingOnly = TRUE)
m <- as.character(args[1])

# Loading the data
dnam <- readRDS("Data/imputed_denoised_450K")
covars <- readRDS("Data/covariates.rds")
print(all(rownames(dnam) == rownames(covars)))

# Model preparation
xdata <- cbind(dnam, age = covars$age.recr, sex = covars$gender)
penalty <- c(rep(1, ncol(dnam)), 0, 0)
if (m == 2) {
  xdata <- cbind(xdata, smoking = covars$csi)
  penalty <- c(penalty, 0)
}

# Variable selection
set.seed(1)
print(system.time({
  mymodel <- cv.glmnet(
    x = xdata, y = covars$case_lc,
    family = "binomial", alpha = 1, type.measure = "class",
    penalty.factor = penalty
  )
}))

# Saving outputs
saveRDS(mymodel, paste0("Results/cv_lasso_m", m, ".rds"))
