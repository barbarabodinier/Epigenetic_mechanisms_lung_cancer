setwd("../Results/")

library(focus)

args <- commandArgs(trailingOnly = TRUE)
print(args)
m <- as.character(args[1])
ttd <- as.numeric(args[2])
print(ttd)

stab <- readRDS(paste0("Iterations/stability_80_m", m, "_", ttd, "_", 1, ".rds"))

for (i in 2:50) {
  print(i)
  tmp <- readRDS(paste0("Iterations/stability_80_m", m, "_", ttd, "_", i, ".rds"))
  print(system.time({
    stab <- Combine(stab, tmp, include_beta = FALSE)
  }))
}

saveRDS(stab, paste0("stability_80_m", m, "_", ttd, ".rds"))
