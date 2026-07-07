## 02_differential_binding_analysis.R
##
## Differential binding analysis for CUT&RUN peak calls:
##   DiffBind (consensus peaks, counting, normalization)
##     -> edgeR (differential test, via DiffBind)
##     -> ChIPseeker (genomic annotation; mouse OR human block below)
##
## See README.md for package choice rationale.
##
## Usage:
##   Rscript scripts/02_differential_binding_analysis.R

## ------------------------- USER CONFIGURATION -------------------------

OUTDIR <- "/path/to/your/project/CUT_and_RUN/results"

# Sample sheet: one row per sample. Condition/Replicate/bamReads/Peaks
# must line up with the sample names and outputs from
# scripts/01_call_peaks.sh.
samples <- data.frame(
  SampleID = c("treated_rep1", "treated_rep2", "treated_rep3",
               "control_rep1", "control_rep2", "control_rep3"),
  Condition = c("Treated", "Treated", "Treated",
                "Control", "Control", "Control"),
  Replicate = c(1, 2, 3, 1, 2, 3),
  bamReads = file.path(OUTDIR, "05_filtered", paste0(
    c("sample_treated_rep1", "sample_treated_rep2", "sample_treated_rep3",
      "sample_control_rep1", "sample_control_rep2", "sample_control_rep3"),
    "_dedup.bam"
  )),
  Peaks = file.path(OUTDIR, "08_peaks/macs3_cutnrun", paste0(
    c("sample_treated_rep1", "sample_treated_rep2", "sample_treated_rep3",
      "sample_control_rep1", "sample_control_rep2", "sample_control_rep3"),
    "_peaks.narrowPeak"
  )),
  PeakCaller = "narrowPeak",
  stringsAsFactors = FALSE
)

# Which condition is "group1" / "group2" in the contrast, and how results
# are labeled below (gained = higher in group1, lost = higher in group2)
GROUP1_NAME <- "Treated"
GROUP2_NAME <- "Control"

## ------------------------------------------------------------------------

library(DiffBind)
library(ChIPseeker)
library(org.Mm.eg.db)   # swap for org.Hs.eg.db if annotating human data
library(TxDb.Mmusculus.UCSC.mm39.refGene)  # swap for hg38 TxDb below
library(ggplot2)

## --- Build consensus peak set, count reads, normalize ---

dba_obj <- dba(sampleSheet = samples, minOverlap = 1)
dba_obj <- dba.count(dba_obj, bUseSummarizeOverlaps = TRUE)
dba_obj <- dba.normalize(dba_obj)
dba_obj <- dba.contrast(
  dba_obj,
  group1 = dba_obj$masks[[GROUP1_NAME]],
  group2 = dba_obj$masks[[GROUP2_NAME]],
  name1 = GROUP1_NAME,
  name2 = GROUP2_NAME
)
dba_obj <- dba.analyze(dba_obj, method = DBA_EDGER)
dba_obj

## --- Differential results (all tested sites, th = 1) ---
## th = 1 / bUsePval = FALSE returns every tested site rather than only
## significant ones, so you can inspect the full FDR distribution —
## useful for checking near-threshold signal in lower-power experiments.

results <- dba.report(dba_obj, method = DBA_EDGER, contrast = 1,
                       th = 1, bUsePval = FALSE)
results_df <- as.data.frame(results)

cat("Total sites tested:", nrow(results_df), "\n")
cat("Significant FDR<0.05:", sum(results_df$FDR < 0.05), "\n")
cat(sprintf("Gained in %s (FDR<0.05): %d\n", GROUP1_NAME,
            sum(results_df$FDR < 0.05 & results_df$Fold > 0)))
cat(sprintf("Lost in %s (FDR<0.05): %d\n", GROUP1_NAME,
            sum(results_df$FDR < 0.05 & results_df$Fold < 0)))

## Columns to keep in output tables
COLS <- c("seqnames", "start", "end", "SYMBOL", "annotation",
          "Conc", paste0("Conc_", GROUP1_NAME), paste0("Conc_", GROUP2_NAME),
          "Fold", "FDR")

outpath <- file.path(OUTDIR, "diffbind_results")
dir.create(outpath, showWarnings = FALSE)

## ========================================================================
## ANNOTATION — choose ONE of the two blocks below depending on species.
## Both blocks take the same `results` GRanges object and produce the
## same `anno_df` data frame used in the rest of the script.
## ========================================================================

## --- Mouse annotation (mm39) ---
## Requires: TxDb.Mmusculus.UCSC.mm39.refGene, org.Mm.eg.db
## (both loaded above; install via install_r_packages.R)

txdb_mouse <- TxDb.Mmusculus.UCSC.mm39.refGene
anno <- annotatePeak(
  results,
  tssRegion = c(-3000, 3000),
  TxDb = txdb_mouse,
  annoDb = "org.Mm.eg.db"
)
anno_df <- as.data.frame(anno)

## --- Human annotation (hg38) ---
## Requires: TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db
## (install via install_r_packages.R; also load them in place of the
## mouse libraries above if this is the only block you're running)
##
## library(org.Hs.eg.db)
## library(TxDb.Hsapiens.UCSC.hg38.knownGene)
##
## txdb_human <- TxDb.Hsapiens.UCSC.hg38.knownGene
## anno <- annotatePeak(
##   results,
##   tssRegion = c(-3000, 3000),
##   TxDb = txdb_human,
##   annoDb = "org.Hs.eg.db"
## )
## anno_df <- as.data.frame(anno)

## ========================================================================
## Downstream reporting (species-agnostic — uses `anno_df` from whichever
## annotation block was run above)
## ========================================================================

anno_df_fdr <- anno_df[order(anno_df$FDR), ]

top10_gained_fdr <- head(
  anno_df_fdr[anno_df_fdr$Fold > 0 & anno_df_fdr$FDR < 0.05, ], 10
)
top10_lost_fdr <- head(
  anno_df_fdr[anno_df_fdr$Fold < 0 & anno_df_fdr$FDR < 0.05, ], 10
)

cat(sprintf("\n=== TOP 10 GAINED PEAKS FDR<0.05 (%s > %s) ===\n",
            GROUP1_NAME, GROUP2_NAME))
print(top10_gained_fdr[, COLS], row.names = FALSE)

cat(sprintf("\n=== TOP 10 LOST PEAKS FDR<0.05 (%s < %s) ===\n",
            GROUP1_NAME, GROUP2_NAME))
print(top10_lost_fdr[, COLS], row.names = FALSE)

write.csv(anno_df_fdr[, COLS],
          file = file.path(outpath, "all_peaks_FDR_annotated.csv"),
          row.names = FALSE)
write.csv(top10_gained_fdr[, COLS],
          file = file.path(outpath, "top10_gained_FDR.csv"),
          row.names = FALSE)
write.csv(top10_lost_fdr[, COLS],
          file = file.path(outpath, "top10_lost_FDR.csv"),
          row.names = FALSE)

## --- Peaks ranked by FDR, regardless of significance cutoff ---
## Useful when enrichment is weak and no peaks clear FDR < 0.05 — lets
## you inspect the most-promising sites even if none are formally
## significant (e.g., while troubleshooting antibody/pulldown efficiency).

gained_near_fdr <- anno_df[anno_df$Fold > 0, ]
gained_near_fdr <- gained_near_fdr[order(gained_near_fdr$FDR), ]

cat("\n=== GAINED PEAKS RANKED BY FDR (most significant first) ===\n")
print(head(gained_near_fdr[, COLS], 15), row.names = FALSE)

lost_near_fdr <- anno_df[anno_df$Fold < 0, ]
lost_near_fdr <- lost_near_fdr[order(lost_near_fdr$FDR), ]

cat("\n=== LOST PEAKS RANKED BY FDR (most significant first) ===\n")
print(head(lost_near_fdr[, COLS], 15), row.names = FALSE)

write.csv(head(gained_near_fdr[, COLS], 15),
          file = file.path(outpath, "gained_ranked_byFDR.csv"),
          row.names = FALSE)
write.csv(head(lost_near_fdr[, COLS], 15),
          file = file.path(outpath, "lost_ranked_byFDR.csv"),
          row.names = FALSE)
