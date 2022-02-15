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

# Loading the data
dnam <- readRDS("Data/imputed_denoised_450K")
covars <- readRDS("Data/covariates.rds")
print(all(rownames(dnam) == rownames(covars)))
annot <- readRDS("Data/hm450.rds")

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
  xdata = cbind(covars[, c("age.recr", "gender")]),
  ydata = ydata, niter = niter
)

roclogistic1 <- rocfoldslogistic(
  xdata = cbind(
    dnam[, rownames(summary_m1), drop = FALSE],
    covars[, c("age.recr", "gender")]
  ),
  ydata = ydata, niter = niter
)

roclogistic1bis <- rocfoldslogistic(
  xdata = cbind(
    dnam[, rownames(summary_m1), drop = FALSE],
    covars[, c("age.recr", "gender", "csi")]
  ),
  ydata = ydata, niter = niter
)

roclogistic2 <- rocfoldslogistic(
  xdata = covars[, c("age.recr", "gender", "csi"), drop = FALSE],
  ydata = ydata, niter = niter
)

roclogistic3 <- rocfoldslogistic(
  xdata = cbind(
    dnam[, rownames(summary_m2), drop = FALSE],
    covars[, c("age.recr", "gender", "csi"), drop = FALSE]
  ),
  ydata = ydata, niter = niter
)

# Incremental performances
for (model_id in 1:2) {
  print(model_id)
  selprop <- eval(parse(text = paste0("selprop_m", model_id)))
  all_variables <- (names(selprop)[sort.list(selprop, decreasing = TRUE)])[1:60]
  auc_incremental_mean <- auc_incremental_lower <- auc_incremental_upper <- NULL
  for (k in 1:length(all_variables)) {
    print(k)
    selected <- all_variables[1:k]
    if (model_id == 1) {
      tmproc <- rocfoldslogistic(
        xdata = cbind(
          dnam[, selected, drop = FALSE],
          covars[, c("age.recr", "gender")]
        ),
        ydata = ydata, niter = niter
      )
    } else {
      tmproc <- rocfoldslogistic(
        xdata = cbind(
          dnam[, selected, drop = FALSE],
          covars[, c("age.recr", "gender", "csi")]
        ),
        ydata = ydata, niter = niter
      )
    }
    auc_incremental_mean <- c(auc_incremental_mean, mean(tmproc$AUC))
    auc_incremental_lower <- c(auc_incremental_lower, quantile(tmproc$AUC, probs = q_lower))
    auc_incremental_upper <- c(auc_incremental_upper, quantile(tmproc$AUC, probs = q_upper))
  }
  names(auc_incremental_mean) <- names(auc_incremental_upper) <- names(auc_incremental_lower) <- all_variables

  assign(paste0("auc_incremental_mean_m", model_id), auc_incremental_mean)
  assign(paste0("auc_incremental_lower_m", model_id), auc_incremental_lower)
  assign(paste0("auc_incremental_upper_m", model_id), auc_incremental_upper)
}

# Figure
{
  pdf("Figures/ROC_analyses.pdf", width = 7, height = 13)
  layout(mat = cbind(c(1, 2, 3)), heights = c(2, 1, 1))

  # ROC curves
  par(mar = c(4.5, 5, 1, 1))
  mymodels <- paste0("roclogistic", 0:3)
  quantiles <- c(q_lower, q_upper)
  colors <- c("darkgreen", "navy", "orange", "darkred")
  tcolors <- adjustcolor(colors, alpha.f = 0.1)

  PlotROCCurves(mymodels, niter = niter)
  par(xpd = TRUE)
  par(xpd = FALSE)

  legends <- c(
    "age + sex",
    paste0("age + sex + CpGs (N=29)"),
    "age + sex + CSI",
    paste0("age + sex + CSI + CpGs (N=50)")
  )
  aucs <- NULL
  for (i in 0:3) {
    roclogistic <- eval(parse(text = paste0("roclogistic", i)))
    aucs <- c(aucs, paste0(
      "- AUC = ", formatC(mean(roclogistic$AUC), digits = 2, format = "f"),
      " [", paste0(formatC(quantile(roclogistic$AUC, probs = quantiles), digits = 2, format = "f"), collapse = "-"), "]"
    ))
  }
  mtext(text = "A", at = 1, side = 2, line = 2.5, las = 1, cex = 2.5)
  legend("bottomright",
    lty = 1, lwd = 2, col = colors, cex = 1.3,
    legend = paste(legends, aucs), bty = "n"
  )

  colours <- c("darkred", "dodgerblue4")

  # Base model
  ids <- names(auc_incremental_mean_m1)
  sm_membership <- ifelse(ids %in% sm1, yes = 1, no = 2)
  mycolours <- colours[sm_membership]

  par(mar = c(11, 5, 1, 1))
  plotCI(auc_incremental_mean_m1,
    ui = auc_incremental_upper_m1,
    li = auc_incremental_lower_m1,
    col = mycolours,
    ylim = c(0.5, 1), pch = 18,
    xlab = "", ylab = "AUC", xaxt = "n",
    las = 1, cex.lab = 1.5, sfrac = 0.003
  )
  abline(h = mean(roclogistic0$AUC), col = "forestgreen", lty = 2)
  abline(v = sum(selected_m1) + 0.5, col = "black", lty = 2)
  for (i in 1:length(auc_incremental_mean_m1)) {
    axis(
      side = 1, at = i, las = 2,
      col = mycolours[i], col.axis = mycolours[i], cex.axis = 0.7,
      labels = paste0("+ ", annot[names(auc_incremental_mean_m1), "alt.name"], "-", names(auc_incremental_mean_m1))[i]
    )
    abline(v = i, col = mycolours[i], lty = 3)
  }
  text(x = 60, y = mean(roclogistic0$AUC), labels = "age + sex", adj = c(1, -1), col = "forestgreen")
  mtext(text = "B", at = 1, side = 2, line = 2.5, las = 1, cex = 2.5)

  # Model adjusted on smoking
  ids <- names(auc_incremental_mean_m2)
  sm_membership <- ifelse(ids %in% sm1, yes = 1, no = 2)
  mycolours <- colours[sm_membership]

  par(mar = c(11, 5, 1, 1))
  plotCI(auc_incremental_mean_m2,
    ui = auc_incremental_upper_m2,
    li = auc_incremental_lower_m2,
    col = mycolours,
    ylim = c(0.5, 1), pch = 18,
    xlab = "", ylab = "AUC", xaxt = "n",
    las = 1, cex.lab = 1.5, sfrac = 0.003
  )
  abline(h = mean(roclogistic2$AUC), col = "orange", lty = 2)
  abline(v = sum(selected_m2) + 0.5, col = "black", lty = 4)
  for (i in 1:length(auc_incremental_mean_m2)) {
    axis(
      side = 1, at = i, las = 2,
      col = mycolours[i], col.axis = mycolours[i], cex.axis = 0.7,
      labels = paste0("+ ", annot[names(auc_incremental_mean_m2), "alt.name"], "-", names(auc_incremental_mean_m2))[i]
    )
    abline(v = i, col = mycolours[i], lty = 3)
  }
  text(x = 60, y = mean(roclogistic2$AUC), labels = "age + sex + CSI", adj = c(1, -1), col = "orange")
  mtext(text = "C", at = 1, side = 2, line = 2.5, las = 1, cex = 2.5)
  dev.off()
}
