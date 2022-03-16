rm(list = ls())

library(glmnet)
library(sharp)
library(openxlsx)

# Loading annotation
annot <- readRDS("Data/hm450.rds")

# Loading lists of smoking-related CpGs
sm1 <- readRDS("Data/List of CpGs (S.London MA)/smk_g1")

# Stability selection LASSO results
stab_m1 <- readRDS("Results/stability_80_m1_0.rds")
class(stab_m1) <- "variable_selection"
selected_m1 <- SelectedVariables(stab_m1)
selprop_m1 <- SelectionProportions(stab_m1)
saveRDS(selprop_m1, "Results/selection_proportions_stability_m1.rds")

stab_m2 <- readRDS("Results/stability_80_m2_0.rds")
class(stab_m2) <- "variable_selection"
selected_m2 <- SelectedVariables(stab_m2)
selprop_m2 <- SelectionProportions(stab_m2)
saveRDS(selprop_m2, "Results/selection_proportions_stability_m2.rds")

# Beta coefficients
tmp <- readRDS("Results/stability_beta_80_m1_0.rds")
beta_m1 <- apply(tmp, 2, FUN = function(x) {
  mean(x[x != 0])
})

tmp <- readRDS("Results/stability_beta_80_m2_0.rds")
beta_m2 <- apply(tmp, 2, FUN = function(x) {
  mean(x[x != 0])
})

# Summary
mat <- cbind(sort(selprop_m1[selected_m1 == 1], decreasing = TRUE))
mat <- cbind(annot[rownames(mat), "alt.name"], mat)
mat <- cbind(mat, formatC(beta_m1[rownames(mat)], format = "f", digits = 2))
colnames(mat) <- c("alt.name", "prop", "beta_mean")
write.table(mat,
  row.names = TRUE, col.names = TRUE, quote = FALSE,
  "Results/Summary_stability_m1.txt"
)

mat <- cbind(sort(selprop_m2[selected_m2 == 1], decreasing = TRUE))
mat <- cbind(annot[rownames(mat), "alt.name"], mat)
mat <- cbind(mat, formatC(beta_m2[rownames(mat)], format = "f", digits = 2))
colnames(mat) <- c("alt.name", "prop", "beta_mean")
write.table(mat,
  row.names = TRUE, col.names = TRUE, quote = FALSE,
  "Results/Summary_stability_m2.txt"
)

# Supplementary Tables
mat <- cbind(sort(selprop_m1[selected_m1 == 1], decreasing = TRUE))
mat <- cbind(
  rownames(mat),
  formatC(beta_m1[rownames(mat)], format = "f", digits = 2),
  mat,
  paste0("chr", annot[rownames(mat), c("chr")]),
  annot[rownames(mat), c("alt.name")],
  ifelse(rownames(mat) %in% sm1, yes = "related", no = "unrelated")
)
colnames(mat) <- c(
  "CpG site", "Average beta", "Selection proportion",
  "Chromosome", "Gene", "Smoking"
)
write.xlsx(as.data.frame(mat),
  "Results/Supplementary_table1.xlsx",
  overwrite = TRUE
)

mat <- cbind(sort(selprop_m2[selected_m2 == 1], decreasing = TRUE))
mat <- cbind(
  rownames(mat),
  formatC(beta_m1[rownames(mat)], format = "f", digits = 2),
  mat,
  paste0("chr", annot[rownames(mat), c("chr")]),
  annot[rownames(mat), c("alt.name")],
  ifelse(rownames(mat) %in% sm1, yes = "related", no = "unrelated")
)
colnames(mat) <- c(
  "CpG site", "Average beta", "Selection proportion",
  "Chromosome", "Gene", "Smoking"
)
write.xlsx(as.data.frame(mat),
  "Results/Supplementary_table2.xlsx",
  overwrite = TRUE
)
