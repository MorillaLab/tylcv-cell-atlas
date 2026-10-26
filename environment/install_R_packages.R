# R packages for the optional limma cells of notebook 02 (11b, 11b-2)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR"), update = FALSE, ask = FALSE)
