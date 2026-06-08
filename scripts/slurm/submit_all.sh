#!/bin/bash
# -----------------------------------------------------------------------------
# Submit the full alignment -> sort/index -> SV-calling pipeline as a chain
# of SLURM jobs with --dependency=afterok so each step waits for the previous.
#
# Usage:
#   bash submit_all.sh <sample_name> <reads.fastq[.gz]> <reference.fa>
# -----------------------------------------------------------------------------

set -euo pipefail

SAMPLE="${1:?sample name required}"
READS="${2:?reads fastq required}"
REF="${3:?reference fasta required}"

mkdir -p logs "results/${SAMPLE}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Submitting NGMLR alignment for ${SAMPLE}"
JID1=$(sbatch --parsable "${SCRIPT_DIR}/01_align_ngmlr.sh" "${SAMPLE}" "${READS}" "${REF}")
echo "  jobid: ${JID1}"

echo "Submitting sort+index (depends on ${JID1})"
JID2=$(sbatch --parsable --dependency=afterok:"${JID1}" \
    "${SCRIPT_DIR}/02_sort_index.sh" "${SAMPLE}")
echo "  jobid: ${JID2}"

echo "Submitting Sniffles2 (depends on ${JID2})"
JID3=$(sbatch --parsable --dependency=afterok:"${JID2}" \
    "${SCRIPT_DIR}/03_sniffles.sh" "${SAMPLE}" "${REF}")
echo "  jobid: ${JID3}"

echo
echo "Pipeline submitted. Track with:"
echo "  squeue -u \$USER"
echo "  tail -f logs/*_${JID1}.out logs/*_${JID2}.out logs/*_${JID3}.out"
