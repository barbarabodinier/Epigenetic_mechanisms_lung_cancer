rm(list = ls())

library(glmnet)
library(sharp)
library(pROC)
library(openxlsx)
library(plotrix)

# Loading annotations
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

# Beta coefficients
tmp <- readRDS("Results/stability_beta_80_m1_0.rds")
beta_m1 <- apply(tmp, 2, FUN = function(x) {
  mean(x[x != 0])
})

tmp <- readRDS("Results/stability_beta_80_m2_0.rds")
beta_m2 <- apply(tmp, 2, FUN = function(x) {
  mean(x[x != 0])
})

# Visualisation of selection proportions
colours <- c("darkred", "dodgerblue4")
for (model_id in c("m1", "m2")) {
  selected <- eval(parse(text = paste0("selected_", model_id)))
  ids <- names(which(selected == 1))
  sm_membership <- ifelse(ids %in% sm1, yes = 1, no = 2)
  mycolours <- colours[sm_membership]
  names(mycolours) <- ids
  assign(paste0("mycolours_", model_id), mycolours)
}

{
  pdf("Figures/Selection_proportions.pdf", width = 15, height = 7.5)
  par(mar = c(12, 5, 1, 1), mfrow = c(1, 2))
  myletters <- LETTERS[1:2]
  names(myletters) <- c("m1", "m2")
  for (model_id in c("m1", "m2")) {
    selected <- eval(parse(text = paste0("selected_", model_id)))
    selprop <- eval(parse(text = paste0("selprop_", model_id)))
    beta <- eval(parse(text = paste0("beta_", model_id)))
    mycolours <- eval(parse(text = paste0("mycolours_", model_id)))
    stab <- eval(parse(text = paste0("stab_", model_id)))
    ordered_selprop <- sort(selprop[which(selected == 1)], decreasing = TRUE)
    plot(seq(1, sum(selected == 1)),
      ordered_selprop,
      type = "h", lwd = 7, lend = 1,
      ylim = c(0, 1.15), xaxt = "n", yaxt = "n", bty = "n",
      xlab = "", ylab = "Selection proportion", cex.lab = 1.5,
      panel.first = c(
        abline(h = seq(0, 1, by = 0.1), lty = 3, col = "grey"),
        abline(h = Argmax(stab)[2], lty = 2, col = "black")
      ),
      col = mycolours[names(ordered_selprop)]
    )
    for (i in 1:length(ordered_selprop)) {
      text(i, ordered_selprop[i] + 0.01,
        srt = 90, cex = 0.8, adj = 0,
        col = mycolours[names(ordered_selprop)[i]],
        labels = eval(parse(text = paste0(
          "expression(beta~'='~'",
          formatC(beta[names(ordered_selprop)[i]], format = "f", digits = 2), "')"
        )))
      )
    }
    mtext(side = 2, at = 1, text = myletters[model_id], las = 1, line = 2.7, cex = 3)
    axis(side = 2, at = seq(0, 1, by = 0.1), las = 1)
    ids <- names(ordered_selprop)
    mylegend <- paste0(annot[ids, "alt.name"], "-", ids)
    axis(side = 1, at = c(1, length(ordered_selprop)), labels = NA)
    for (i in 1:length(ordered_selprop)) {
      axis(
        side = 1, at = i,
        col.axis = mycolours[names(ordered_selprop)[i]],
        labels = mylegend[i], las = 2, cex.axis = 0.8
      )
    }
  }
  dev.off()
}

# Comparison of selection proportions
sum(names(selected_m2)[selected_m2 == 1] %in% names(selected_m1)[selected_m1 == 1])
ids <- sort(unique(c(names(selected_m1)[selected_m1 == 1], names(selected_m2)[selected_m2 == 1])))
ids <- c(ids[ids %in% sm1], ids[!ids %in% sm1])
ids_overlap <- intersect(names(selected_m2)[selected_m2 == 1], names(selected_m1)[selected_m1 == 1])

sm_membership <- ifelse(ids %in% sm1, yes = 1, no = 2)
mycolours <- colours[sm_membership]
pchs <- c(15, 19, 19)
mypchs <- pchs[sm_membership]

{
  pdf("Figures/Scatter_plot.pdf", width = 15, height = 7.5)
  par(mar = c(5, 5, 1, 1), mfrow = c(1, 2))
  plot(selprop_m1[ids], selprop_m2[ids],
    col = mycolours,
    pch = mypchs,
    cex = 0.5,
    xlim = c(0, 1), ylim = c(0, 1),
    xlab = "Selection proportion (base model)",
    ylab = "Selection proportion (model adjusted on CSI)",
    cex.lab = 1.5,
    las = 1,
    panel.first = c(
      abline(v = Argmax(stab_m1)[2], lty = 2),
      abline(h = Argmax(stab_m2)[2], lty = 2)
    )
  )
  mtext(side = 2, at = 1, text = "C", las = 1, line = 2.7, cex = 3)
  text(selprop_m1[ids], selprop_m2[ids],
    offset = 0.3,
    labels = 1:length(ids), pos = 3, col = mycolours
  )
  plot(beta_m1[ids], beta_m2[ids],
    col = mycolours,
    pch = mypchs,
    cex = 0.5,
    xlab = expression("Average" ~ beta ~ "(base model)"),
    ylab = expression("Average" ~ beta ~ "(model adjusted on CSI)"),
    cex.lab = 1.5,
    las = 1,
    panel.first = c(
      abline(v = 0, lty = 2),
      abline(h = 0, lty = 2),
      abline(0, 1, lty = 2)
    )
  )
  mtext(side = 2, at = max(beta_m2[ids]), text = "D", las = 1, line = 2.8, cex = 3)
  text(beta_m1[ids], beta_m2[ids],
    offset = 0.3,
    labels = 1:length(ids), pos = 3, col = mycolours
  )
  dev.off()
}

plotname <- "Figures/Scatter_plot_legend.pdf"
mylegend <- paste0(1:length(ids), ": ", annot[ids, "alt.name"], "-", ids)
{
  pdf(plotname, width = 18, height = 6)
  plot.new()
  legend("top",
    bty = "n",
    legend = sapply(mylegend, FUN = function(x) {
      if (gsub(".*-", "", x) %in% ids_overlap) {
        eval(parse(text = paste0("expression(bold('", x, "'))")))
      } else {
        eval(parse(text = paste0("expression('", x, "')")))
      }
    }),
    text.col = mycolours, ncol = 5, y.intersp = 1.5
  )
  dev.off()
}
system(paste("pdfcrop --margin 10", plotname, plotname))
