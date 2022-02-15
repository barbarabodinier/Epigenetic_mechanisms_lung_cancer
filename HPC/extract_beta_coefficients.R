setwd("../Results/")

library(focus)

args <- commandArgs(trailingOnly = TRUE)
m <- as.character(args[1])
ttd <- as.numeric(args[2])

stab <- readRDS(paste0("stability_80_m", m, "_", ttd, ".rds"))
argmax_id <- ArgmaxId(stab)[1]
print(argmax_id)

Beta <- matrix(NA, nrow = 1000, ncol = ncol(stab$selprop))
colnames(Beta) <- colnames(stab$selprop)
print(dim(Beta))
id <- 0
for (i in 1:50) {
  print(i)
  tmp <- readRDS(paste0("Iterations/stability_80_m", m, "_", ttd, "_", i, ".rds"))
  print(dim(tmp$Beta))
  Beta[(id + 1):(id + tmp$params$K), ] <- t(tmp$Beta[argmax_id, colnames(Beta), ])
  id <- id + tmp$params$K
}

saveRDS(Beta, paste0("stability_beta_80_m", m, "_", ttd, ".rds"))
