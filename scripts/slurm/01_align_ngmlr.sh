#!/bin/bash
#SBATCH --job-name=ngmlr_align
#SBATCH --output=logs/ngmlr_%j.out
#SBATCH --error=logs/ngmlr_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=long
# -----------------------------------------------------------------------------
# Align Oxford Nanopore reads to a reference genome with NGMLR.
#
# Usage (standalone):
#   sbatch 01_align_ngmlr.sh <sample_name> <reads.fastq[.gz]> <reference.fa>
#
# Or driven by submit_all.sh (recommended).
# -----------------------------------------------------------------------------

set -euo pipefail

SAMPLE="${1:?sample name required}"
READS="${2:?reads fastq required}"
REF="${3:?reference fasta required}"

OUTDIR="${OUTDIR:-results/${SAMPLE}}"
mkdir -p "${OUTDIR}" logs

THREADS="${SLURM_CPUS_PER_TASK:-8}"

echo "[$(date)] NGMLR alignment | sample=${SAMPLE} | threads=${THREADS}"

ngmlr \
    -t "${THREADS}" \
    -r "${REF}" \
    -q "${READS}" \
    -o "${OUTDIR}/${SAMPLE}.sam" \
    -x ont

echo "[$(date)] Converting SAM -> BAM"
samtools view -@ "${THREADS}" -bS "${OUTDIR}/${SAMPLE}.sam" \
    > "${OUTDIR}/${SAMPLE}.unsorted.bam"

rm "${OUTDIR}/${SAMPLE}.sam"

echo "[$(date)] Done. Unsorted BAM: ${OUTDIR}/${SAMPLE}.unsorted.bam"
