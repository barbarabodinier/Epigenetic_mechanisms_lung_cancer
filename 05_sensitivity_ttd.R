rm(list = ls())

library(glmnet)
library(focus)

# Loading annotation
annot <- readRDS("Data/hm450.rds")

ttd <- 5

# Stability selection LASSO results
stab_m1 <- readRDS("Results/stability_80_m1_0.rds")
class(stab_m1) <- "variable_selection"
selected_m1 <- SelectedVariables(stab_m1)
selprop_m1 <- SelectionProportions(stab_m1)
selprop_m1_0 <- selprop_m1
selected_m1_0 <- selected_m1
pi_m1_0 <- Argmax(stab_m1)[2]

stab_m2 <- readRDS("Results/stability_80_m2_0.rds")
class(stab_m2) <- "variable_selection"
selected_m2 <- SelectedVariables(stab_m2)
selprop_m2 <- SelectionProportions(stab_m2)
selprop_m2_0 <- selprop_m2
selected_m2_0 <- selected_m2
pi_m2_0 <- Argmax(stab_m2)[2]

stab_m1 <- readRDS(paste0("Results/stability_80_m1_", ttd, ".rds"))
class(stab_m1) <- "variable_selection"
selected_m1 <- SelectedVariables(stab_m1)
selprop_m1 <- SelectionProportions(stab_m1)
selprop_m1_1 <- selprop_m1
selected_m1_1 <- selected_m1
pi_m1_1 <- Argmax(stab_m1)[2]

stab_m2 <- readRDS(paste0("Results/stability_80_m2_", ttd, ".rds"))
class(stab_m2) <- "variable_selection"
selected_m2 <- SelectedVariables(stab_m2)
selprop_m2 <- SelectionProportions(stab_m2)
selprop_m2_1 <- selprop_m2
selected_m2_1 <- selected_m2
pi_m2_1 <- Argmax(stab_m2)[2]

{
  png(paste0("Figures/Scatter_plot_ttd_", ttd, "_year.png"),
    width = 1200, height = 600
  )
  par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
  plot(selprop_m1_0, selprop_m1_1,
    xlim = c(0, 1), ylim = c(0, 1),
    pch = 19, col = "navy", las = 1, cex = 0.5,
    cex.lab = 1.5, cex.main = 2,
    panel.first = c(
      abline(0, 1, lty = 2),
      abline(v = pi_m1_0, lty = 3, col = "red"),
      abline(h = pi_m1_1, lty = 3, col = "red")
    ),
    main = "Base model",
    xlab = "Full sample",
    ylab = paste0("Time to diagnosis above ", ttd, " year")
  )
  plot(selprop_m2_0, selprop_m2_1,
    xlim = c(0, 1), ylim = c(0, 1),
    pch = 19, col = "navy", las = 1, cex = 0.5,
    cex.lab = 1.5, cex.main = 2,
    panel.first = c(
      abline(0, 1, lty = 2),
      abline(v = pi_m2_0, lty = 3, col = "red"),
      abline(h = pi_m2_1, lty = 3, col = "red")
    ),
    main = "CSI-adjusted model",
    xlab = "Full sample",
    ylab = paste0("Time to diagnosis above ", ttd, " year")
  )
  dev.off()
}

table(selected_m1_0, selected_m1_1)
table(selected_m2_0, selected_m2_1)
