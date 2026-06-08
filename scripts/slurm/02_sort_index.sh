#!/bin/bash
#SBATCH --job-name=sort_index
#SBATCH --output=logs/sort_index_%j.out
#SBATCH --error=logs/sort_index_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=short
# -----------------------------------------------------------------------------
# Sort and index an unsorted BAM with samtools.
#
# Usage:
#   sbatch 02_sort_index.sh <sample_name>
# Expects:
#   results/<sample>/<sample>.unsorted.bam
# Produces:
#   results/<sample>/<sample>.sorted.bam
#   results/<sample>/<sample>.sorted.bam.bai
# -----------------------------------------------------------------------------

set -euo pipefail

SAMPLE="${1:?sample name required}"
OUTDIR="${OUTDIR:-results/${SAMPLE}}"

IN_BAM="${OUTDIR}/${SAMPLE}.unsorted.bam"
SORTED_BAM="${OUTDIR}/${SAMPLE}.sorted.bam"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

echo "[$(date)] Sorting ${IN_BAM}"
samtools sort -@ "${THREADS}" -m 3G -o "${SORTED_BAM}" "${IN_BAM}"

echo "[$(date)] Indexing ${SORTED_BAM}"
samtools index -@ "${THREADS}" "${SORTED_BAM}"

# Quick QC summary
echo "[$(date)] flagstat:"
samtools flagstat "${SORTED_BAM}" | tee "${OUTDIR}/${SAMPLE}.flagstat.txt"

# Remove the unsorted file to save space (comment out to keep)
rm -f "${IN_BAM}"

echo "[$(date)] Done. Sorted BAM ready for Sniffles2."
