# CUT&RUN Differential Binding Analysis Pipeline

A reproducible pipeline for calling peaks and identifying differentially
bound regions from CUT&RUN sequencing data, with downstream annotation
against a reference genome. Built and tested on a transcription factor
(ESR1) CUT&RUN experiment comparing hormone-treated vs. vehicle-control
conditions in mouse, but written to generalize to any two-condition
CUT&RUN comparison in mouse or human.

## Pipeline overview

```
Raw reads
   │  (FastQC, adapter trimming, alignment, deduplication — upstream of this repo)
   ▼
Deduplicated, filtered BAM files (target + IgG control)
   │
   ▼
01_call_peaks.sh          → MACS3 peak calling (per-sample narrowPeak files)
   │
   ▼
02_differential_binding_analysis.R
   ├── DiffBind: consensus peak set, read counting, normalization
   ├── edgeR (via DiffBind): differential binding contrast
   └── ChIPseeker: genomic annotation (mouse block + human block)
   ▼
Annotated differential binding tables (CSV) + summary stats
```

This repo picks up **after** raw read QC/trimming/alignment/deduplication
(e.g., FastQC → Trimmomatic → Bowtie2 → Picard `MarkDuplicates`), since that
part of the workflow is largely identical across CUT&RUN/ChIP-style
experiments and is well covered by existing pipelines (e.g., nf-core/cutandrun).
The scripts here assume you're starting from deduplicated, filtered BAM
files for each target sample and a matched IgG (or other isotype) control.

## Repository structure

```
cutnrun-pipeline/
├── README.md
├── LICENSE
├── environment.yml              # conda environment for MACS3 + supporting CLI tools
├── install_r_packages.R         # Bioconductor/CRAN package installation for the R step
├── .gitignore
└── scripts/
    ├── 01_call_peaks.sh          # MACS3 peak calling
    └── 02_differential_binding_analysis.R   # DiffBind + edgeR + ChIPseeker (mouse & human blocks)
```

Both scripts keep all user-specific settings (paths, sample names, condition
labels, genome build) in a single `USER CONFIGURATION` block near the top of
the file — there's no separate config file to keep in sync, and no extra
YAML-parsing dependency. Edit that block directly for your own data.

## Requirements

- **Command-line tools:** MACS3, samtools (see `environment.yml`)
- **R (≥4.2):** DiffBind, ChIPseeker, edgeR, ggplot2, plus a genome-specific
  `TxDb` and `org.db` pair (see `install_r_packages.R`)

```bash
conda env create -f environment.yml
conda activate cutnrun-pipeline
Rscript install_r_packages.R
```

## Usage

1. Edit the `USER CONFIGURATION` block at the top of
   `scripts/01_call_peaks.sh` (sample names, BAM directory, IgG control
   sample, genome size code) to point at your own data.
2. Call peaks:
   ```bash
   bash scripts/01_call_peaks.sh
   ```
3. Edit the `USER CONFIGURATION` block at the top of
   `scripts/02_differential_binding_analysis.R` (sample sheet, output
   directory), then run the differential binding analysis, choosing the
   mouse or human annotation block near the bottom of the script:
   ```bash
   Rscript scripts/02_differential_binding_analysis.R
   ```

Outputs are written to `<OUTDIR>/08_peaks/` (peak calls) and
`<OUTDIR>/diffbind_results/` (annotated differential binding tables).

## Why these tools? (rationale for package choices)

**MACS3 for peak calling, rather than SEACR.**
SEACR was purpose-built for CUT&RUN/CUT&Tag data and is a reasonable
default when no input/IgG control is available, since it calls peaks from
signal thresholds without needing a background model. However, this
experiment *does* have a matched IgG control per replicate, and the
downstream step is a **quantitative differential binding comparison**
(DiffBind/edgeR) rather than a simple presence/absence peak call. MACS3's
model-based approach (`-c` control subtraction, BAMPE mode for paired-end
fragments) gives peak sets with signal/enrichment statistics that behave
well as input to count-based differential frameworks, and it's the peak
caller DiffBind's documentation and most published CUT&RUN differential
workflows are built around. `--nomodel --extsize 200` is used because
CUT&RUN fragments don't have the same shift/extension profile that MACS's
default ChIP-seq fragment size model assumes, so the fragment size is set
explicitly rather than estimated.

**DiffBind for consensus peak building, counting, and normalization.**
DiffBind is the standard Bioconductor package for turning a set of
per-sample peak calls and BAMs into a single consensus binding matrix,
handling read counting (`dba.count`) and normalization (`dba.normalize`)
in a way that's consistent with how the differential test is run
downstream. It also keeps the whole comparison (peaks, counts,
normalization, contrast, and test) reproducible from one object.

**edgeR (via `DBA_EDGER`) for the differential test.**
DiffBind supports both edgeR and DESeq2 as its underlying test engines.
edgeR's negative-binomial model with empirical Bayes dispersion shrinkage
is a good fit here given the small number of replicates per group (n=3),
where per-gene/per-peak dispersion estimates are otherwise noisy. DESeq2
is an equally defensible alternative and swapping `method = DBA_EDGER` for
`method = DBA_DESEQ2` in `dba.analyze()` is a one-line change if you'd
rather cross-check with both.

**ChIPseeker + TxDb/org.db pair for annotation.**
ChIPseeker annotates each peak/region against a transcript database
(nearest gene, distance to TSS, genomic feature — promoter, exon, intron,
intergenic, etc.) and maps IDs to gene symbols via an organism annotation
package. This repo uses:
- Mouse: `TxDb.Mmusculus.UCSC.mm39.refGene` + `org.Mm.eg.db`
- Human: `TxDb.Hsapiens.UCSC.hg38.knownGene` + `org.Hs.eg.db`

Swap in the appropriate `TxDb`/`org.db` pair (and matching genome build)
for other species or genome builds — see the "Human annotation" block in
`02_differential_binding_analysis.R` for a working example.

## Notes / caveats

- Peak counts and enrichment are highly dependent on antibody/pulldown
  efficiency. If one condition shows dramatically higher peak counts
  than expected (e.g., an order of magnitude difference between
  replicates of the same condition), check enrichment/QC metrics (FRiP,
  IgG background) before interpreting differential results — this can
  indicate an antibody or pulldown issue rather than a biological signal.
- `th = 1, bUsePval = FALSE` in `dba.report()` returns **all** tested
  sites (not just significant ones) so you can inspect the full FDR
  distribution — useful when checking for near-threshold signal in
  lower-power experiments. Filter to your significance threshold
  (e.g., FDR < 0.05) before treating results as a final gene list.

## License

MIT — see `LICENSE`.
