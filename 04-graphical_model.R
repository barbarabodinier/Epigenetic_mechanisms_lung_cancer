rm(list = ls())

library(openxlsx)
library(RColorBrewer)
library(sharp)
library(igraph)

# Loading the data
dnam <- readRDS("Data/imputed_denoised_450K")
covars <- readRDS("Data/covariates.rds")
sm1 <- readRDS("Data/List of CpGs (S.London MA)/smk_g1")

# Looping over sheet IDs for CpGs selected in univariate/LASSO analyses
for (mysheet in 2) {
  print(mysheet)

  # Loading the sets of CpGs
  if (mysheet == 1) {
    m1_bonf <- read.xlsx("Data/For network (Nov. 2020).xlsx", sheet = mysheet)
    m1_bonf <- rbind(colnames(m1_bonf), m1_bonf)
    m1_bonf <- m1_bonf[, 1]
    cpg <- dnam[, m1_bonf]
  }
  if (mysheet == 2) {
    m1_bonf <- read.table("Results/Summary_stability_m1.txt")
    m1_bonf <- rownames(m1_bonf)
    cpg <- dnam[, m1_bonf]
  }
  if (mysheet == 3) {
    m1_bonf <- read.table("Results/Selected_stability_m2.txt")
    m1_bonf <- m1_bonf[, 1]
    cpg <- dnam[, m1_bonf]
  }

  # Loading annotation file
  cpgannot <- readRDS("Data/hm450.rds")
  cpgannot <- cpgannot[m1_bonf, ]
  print(all(rownames(cpg) == rownames(covars)))
  x <- cbind(cpg, CSI = covars$csi)
  print(dim(x))
  x <- x[covars$case_lc == 1, ]


  ### Network estimation

  # PFER_thr=10
  PFER_thr <- 10
  out <- GraphicalModel(xdata = x, PFER_thr = PFER_thr, refine_calib_grid = FALSE)
  saveRDS(out, paste0("Results/calibration_", mysheet, "_", ncol(cpg), "_cpg_csi_in_PFER", PFER_thr, ".rds"))
  out <- readRDS(paste0("Results/calibration_", mysheet, "_", ncol(cpg), "_cpg_csi_in_PFER", PFER_thr, ".rds"))

  # Calibration plot
  {
    pdf(paste0("Figures/calibration_heatmap_", mysheet, "_", ncol(cpg), "_cpg_csi_in_PFER", PFER_thr, ".pdf"))
    par(mar = c(5, 5, 5, 5))
    CalibrationPlot(out)
    dev.off()
  }

  adjacency <- Adjacency(out)


  ### Network visualisation

  # Colour by chromosome
  colors <- c(brewer.pal(12, name = "Set3"), brewer.pal(n = 8, name = "Set2"), brewer.pal(n = 4, name = "Set1"), "darkred")
  chr_number <- as.character(cpgannot$chr)
  chr_number[chr_number == "X"] <- 23
  chr_number[chr_number == "Y"] <- 24
  chr_number <- as.numeric(chr_number)
  chr_number <- c(chr_number, 25)

  # Labels with gene name
  colnames(adjacency) <- rownames(adjacency) <- c(rownames(cpgannot), "CSI")
  node_label <- c(paste(cpgannot$alt.name, "\n", rownames(cpgannot)), "CSI")

  # Making igraph object
  gg <- Graph(adjacency = adjacency, node_colour = colors[chr_number], node_label = node_label)
  V(gg)$size <- sqrt(3 * V(gg)$size)
  V(gg)$label.color[V(gg)$name == "CSI"] <- "white"
  V(gg)$shape <- rep("circle", length(V(gg)$name))
  V(gg)$shape[V(gg)$name == "CSI"] <- "square"

  myasp <- 0.7
  {
    pdf(paste0("Figures/network_", mysheet, "_", ncol(cpg), "_cpg_csi_in_PFER", PFER_thr, ".pdf"),
      height = myasp * 12, width = 12
    )
    par(mar = rep(0, 4))
    set.seed(1)
    plot(gg, layout = layout_with_fr(gg), asp = myasp)
    dev.off()
  }

  # Colour by smoking group
  colors <- c("tomato3", "dodgerblue", "black")
  sm_membership <- ifelse(colnames(adjacency) %in% sm1, yes = 1, no = 2)
  sm_membership[length(sm_membership)] <- 3
  nodecolor <- colors[sm_membership]
  names(nodecolor) <- colnames(adjacency)

  # Making igraph object
  gg <- Graph(adjacency = adjacency, node_colour = nodecolor, node_label = node_label, satellites = TRUE)
  V(gg)$size <- sqrt(3 * V(gg)$size)
  V(gg)$label.color[V(gg)$name == "CSI"] <- "white"
  V(gg)$shape <- rep("circle", length(V(gg)$name))
  V(gg)$shape[V(gg)$name == "CSI"] <- "square"

  mycounts <- table(as.factor(nodecolor)[V(gg)$name])

  myasp <- 0.7
  {
    pdf(paste0("Figures/network_", mysheet, "_", ncol(cpg), "_cpg_csi_in_PFER", PFER_thr, "_colour_smoking.pdf"),
      height = myasp * 12, width = 12
    )
    par(mar = rep(0, 4))
    set.seed(1)
    if (mysheet == 1) {
      plot(gg, layout = layout_with_fr(gg), asp = myasp)
    } else {
      plot(gg, layout = layout_with_kk(gg), asp = myasp)
    }
    legend("topleft",
      pch = 19, col = colors, bty = "n",
      legend = c(
        paste(mycounts["tomato3"], "smoking-related CpG sites"),
        # paste(mycounts["darkgoldenrod3"], "Smk-g2-CpG(s)"),
        paste(mycounts["dodgerblue"], "smoking-unrelated CpG sites")
      )
    )
    dev.off()
  }
}
