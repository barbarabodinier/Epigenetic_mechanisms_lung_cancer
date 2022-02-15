rm(list = ls())

library(glmnet)
library(focus)

# library(devtools)
# unzip("Scripts/pROC-master.zip")
# install("Scripts/pROC-master/")
library(pROC)

library(openxlsx)
library(plotrix)
source("Scripts/functions.R")

for (sex in c("Male", "Female")) {
  print(sex)

  # Loading the data
  dnam <- readRDS("Data/imputed_denoised_450K")
  covars <- readRDS("Data/covariates.rds")
  print(all(rownames(dnam) == rownames(covars)))
  annot <- readRDS("Data/hm450.rds")

  # Sex-specific analyses
  ids <- which(covars$gender == sex)
  covars <- covars[ids, ]
  dnam <- dnam[rownames(covars), ]

  # Loading lists of smoking-related CpGs
  sm1 <- readRDS("Data/List of CpGs (S.London MA)/smk_g1")

  # Stability selection LASSO results
  stab_m1 <- readRDS("Results/stability_80_m1_0.rds")
  class(stab_m1) <- "variable_selection"
  selected_m1 <- SelectedVariables(stab_m1)
  selprop_m1 <- SelectionProportions(stab_m1)

  stab_m2 <- readRDS("Results/stability_80_m2_0.rds")
  class(stab_m2) <- "variable_selection"
  selected_m2 <- SelectedVariables(stab_m2)
  selprop_m2 <- SelectionProportions(stab_m2)

  # Loading results summary
  summary_m1 <- read.table("Results/Summary_stability_m1.txt")
  summary_m2 <- read.table("Results/Summary_stability_m2.txt")

  # ROC analyses
  niter <- 200
  q_lower <- 0.05
  q_upper <- 0.95

  ydata <- covars$case_lc
  names(ydata) <- rownames(covars)

  roclogistic0 <- rocfoldslogistic(
    xdata = cbind(covars[, c("age.recr"), drop = FALSE]),
    ydata = ydata, niter = niter
  )

  roclogistic1 <- rocfoldslogistic(
    xdata = cbind(
      dnam[, rownames(summary_m1), drop = FALSE],
      covars[, c("age.recr"), drop = FALSE]
    ),
    ydata = ydata, niter = niter
  )

  roclogistic2 <- rocfoldslogistic(
    xdata = covars[, c("age.recr", "csi"), drop = FALSE],
    ydata = ydata, niter = niter
  )

  roclogistic3 <- rocfoldslogistic(
    xdata = cbind(
      dnam[, rownames(summary_m2), drop = FALSE],
      covars[, c("age.recr", "csi"), drop = FALSE]
    ),
    ydata = ydata, niter = niter
  )


  # ROC curves
  {
    pdf(paste0("Figures/ROC_analyses_", tolower(sex), ".pdf"), width = 7, height = 7)
    par(mar = c(4.5, 5, 1, 1))
    mymodels <- paste0("roclogistic", 0:3)
    quantiles <- c(q_lower, q_upper)
    colors <- c("darkgreen", "navy", "orange", "darkred")
    tcolors <- adjustcolor(colors, alpha.f = 0.1)

    PlotROCCurves(mymodels, niter = niter)
    par(xpd = TRUE)
    par(xpd = FALSE)

    legends <- c(
      "age",
      paste0("age + CpGs (N=29)"),
      "age + CSI",
      paste0("age + CSI + CpGs (N=50)")
    )
    aucs <- NULL
    for (i in 0:3) {
      roclogistic <- eval(parse(text = paste0("roclogistic", i)))
      aucs <- c(aucs, paste0(
        "- AUC = ", formatC(mean(roclogistic$AUC), digits = 2, format = "f"),
        " [", paste0(formatC(quantile(roclogistic$AUC, probs = quantiles), digits = 2, format = "f"), collapse = "-"), "]"
      ))
    }
    legend("bottomright",
      lty = 1, lwd = 2, col = colors, cex = 0.9,
      legend = paste(legends, aucs), bty = "n"
    )
    dev.off()
  }
}
