rm(list = ls())

library(focus)

# Loading the data
dnam <- readRDS("Data/imputed_denoised_450K")
covars <- readRDS("Data/covariates.rds")
stab_m1 <- readRDS("Results/stability_80_m1_0.rds")
class(stab_m1) <- "variable_selection"
stab_m2 <- readRDS("Results/stability_80_m2_0.rds")
class(stab_m2) <- "variable_selection"

# Distribution of WBC proportions
wbc <- colnames(covars)[62:68]
par(mfrow = c(4, 2))
for (i in 1:length(wbc)) {
  plot(density(covars[, wbc[i]]))
}
wbc <- wbc[-6] # removing eosinophils with skewed distribution

# Associations between wbc and lung cancer status
pval <- NULL
for (i in 1:length(wbc)) {
  mymodel <- glm(covars$case_lc ~ covars[, wbc[i]], family = "binomial")
  pval <- c(pval, summary(mymodel)$coefficients[2, 4])
}

for (m_id in c(1, 2)) {
  print(m_id)

  # Calibrated models
  stab <- eval(parse(text = paste0("stab_m", m_id)))
  selected <- SelectedVariables(stab)
  selected <- names(selected)[which(selected == 1)]

  # Not adjusted on WBC
  if (m_id == 1) {
    xdata <- cbind(dnam[, selected], age = covars$age.recr, sex = covars$gender)
  } else {
    xdata <- cbind(dnam[, selected], age = covars$age.recr, sex = covars$gender, smoking = covars$csi)
  }
  mymodel <- Recalibrate(
    xdata = xdata,
    ydata = covars$case_lc,
    family = "binomial"
  )

  # Adjusted on WBC
  if (m_id == 1) {
    xdata <- cbind(dnam[, selected],
      age = covars$age.recr, sex = covars$gender,
      covars[, wbc]
    )
  } else {
    xdata <- cbind(dnam[, selected],
      age = covars$age.recr, sex = covars$gender, smoking = covars$csi,
      covars[, wbc]
    )
  }
  mymodel2 <- Recalibrate(
    xdata = xdata,
    ydata = covars$case_lc,
    family = "binomial"
  )

  {
    pdf(paste0("Figures/Scatter_plot_wbc_m", m_id, ".pdf"))
    par(mar = c(5, 5, 1, 1))
    plot(coefficients(mymodel)[-1],
      coefficients(mymodel2)[names(coefficients(mymodel)[-1])],
      pch = 19, col = "navy", las = 1,
      xlab = ifelse(m_id == 1, yes = "Base model", no = "CSI-adjusted model"),
      ylab = paste0(
        ifelse(m_id == 1, yes = "Base model", no = "CSI-adjusted model"),
        " (adjusted on WBC)"
      ),
      cex.lab = 1.5,
      panel.first = c(
        abline(h = 0, lty = 2),
        abline(v = 0, lty = 2),
        abline(0, 1, lty = 2)
      )
    )
    dev.off()
  }
}
