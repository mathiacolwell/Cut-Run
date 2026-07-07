## install_r_packages.R

## cute way to find all the packges you need and why

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# --- Core differential binding + annotation packages ---
BiocManager::install(c(
  "DiffBind",   # consensus peak building, read counting, normalization, contrasts
  "edgeR",      # differential test engine used via DiffBind (DBA_EDGER)
  "ChIPseeker"  # peak-to-gene annotation, genomic feature classification
), update = FALSE, ask = FALSE)

# --- Mouse annotation resources (mm39) ---
BiocManager::install(c(
  "TxDb.Mmusculus.UCSC.mm39.refGene",
  "org.Mm.eg.db"
), update = FALSE, ask = FALSE)

# --- Human annotation resources (hg38) ---
# Only needed if you're running the human annotation block in
# 02_differential_binding_analysis.R
BiocManager::install(c(
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db"
), update = FALSE, ask = FALSE)

# --- Plotting ---
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

cat("\nAll packages installed. Session info:\n")
print(sessionInfo())
