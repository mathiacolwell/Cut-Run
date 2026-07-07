#!/usr/bin/env bash
##
## Step1_Call_peaks.sh
##
## Calls peaks for each target sample against a shared IgG (or other
## isotype) control using MACS3. See README.md for the rationale behind
## using MACS3 (with an input control) rather than SEACR here.
##
## Usage:
##   bash scripts/01_call_peaks.sh

set -euo pipefail

## ------------------------- USER CONFIGURATION -------------------------

# Root directory containing your deduplicated, filtered BAM files
# (e.g., output of Picard MarkDuplicates + filtering).
BASE_DIR="/path/to/your/project/CUT_and_RUN"

# Sub-paths relative to BASE_DIR
BAM_DIR="${BASE_DIR}/results/05_filtered"
PEAK_OUTDIR="${BASE_DIR}/results/08_peaks/macs3_cutnrun"
LOG_DIR="${BASE_DIR}/results/logs"

# Target sample BAM prefixes (without "_dedup.bam"), one per replicate
SAMPLES=(
  "sample_treated_rep1"
  "sample_treated_rep2"
  "sample_treated_rep3"
  "sample_control_rep1"
  "sample_control_rep2"
  "sample_control_rep3"
)

# Shared IgG / isotype control BAM prefix (without "_dedup.bam")
IGG_CONTROL="sample_IgG"

# Genome size code for MACS3 (-g): "mm" = mouse, "hs" = human,
# or an effective genome size integer for other organisms
GENOME_SIZE="mm"

# Fragment size used for --extsize. CUT&RUN fragments don't follow MACS's
# default ChIP-seq shift/extension model, so this is set explicitly
# rather than estimated by MACS3.
EXTSIZE=200

## ------------------------------------------------------------------------

mkdir -p "${PEAK_OUTDIR}" "${LOG_DIR}"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Calling peaks: ${SAMPLE}"
    macs3 callpeak \
        -t "${BAM_DIR}/${SAMPLE}_dedup.bam" \
        -c "${BAM_DIR}/${IGG_CONTROL}_dedup.bam" \
        -f BAMPE \
        -g "${GENOME_SIZE}" \
        --nomodel \
        --extsize "${EXTSIZE}" \
        --keep-dup all \
        -q 0.05 \
        --outdir "${PEAK_OUTDIR}" \
        -n "${SAMPLE}" \
        2> "${LOG_DIR}/${SAMPLE}_macs3_cutnrun.log"
done

echo "=== Peak counts ==="
for f in "${PEAK_OUTDIR}"/*_peaks.narrowPeak; do
    echo "$(basename "$f"): $(wc -l < "$f") peaks"
done
