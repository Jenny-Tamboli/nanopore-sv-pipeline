#!/bin/bash
#SBATCH --job-name=sniffles
#SBATCH --output=logs/sniffles_%j.out
#SBATCH --error=logs/sniffles_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=short
# -----------------------------------------------------------------------------
# Call structural variants from a sorted+indexed BAM using Sniffles2.
#
# Usage:
#   sbatch 03_sniffles.sh <sample_name> <reference.fa>
# Expects:
#   results/<sample>/<sample>.sorted.bam
# Produces:
#   results/<sample>/<sample>.sniffles.vcf
# -----------------------------------------------------------------------------

set -euo pipefail

SAMPLE="${1:?sample name required}"
REF="${2:?reference fasta required}"
OUTDIR="${OUTDIR:-results/${SAMPLE}}"

SORTED_BAM="${OUTDIR}/${SAMPLE}.sorted.bam"
OUT_VCF="${OUTDIR}/${SAMPLE}.sniffles.vcf"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

echo "[$(date)] Sniffles2 on ${SORTED_BAM}"

sniffles \
    --input "${SORTED_BAM}" \
    --reference "${REF}" \
    --vcf "${OUT_VCF}" \
    --sample-id "${SAMPLE}" \
    --threads "${THREADS}" \
    --minsupport auto \
    --minsvlen 50

echo "[$(date)] Done. VCF: ${OUT_VCF}"

# Optional: immediately run downstream filtering
# python scripts/python/filter_vcf.py "${OUT_VCF}" -o "${OUTDIR}/${SAMPLE}.filtered.vcf"
