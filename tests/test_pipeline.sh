#!/bin/bash
# -----------------------------------------------------------------------------
# End-to-end demo: filter -> CSV -> histograms on the synthetic example VCF.
# Run from the repo root:
#   bash tests/test_pipeline.sh
# -----------------------------------------------------------------------------

set -euo pipefail

SAMPLE="test_sample"
OUTDIR="results/${SAMPLE}"
HISTDIR="${OUTDIR}/histplots"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"

mkdir -p "${OUTDIR}" "${HISTDIR}"

echo "==> 1. Filter VCF (PASS + PRECISE, standard chromosomes only)"
python scripts/python/filter_vcf.py \
    tests/data/example.vcf \
    -o "${OUTDIR}/${SAMPLE}.filtered.vcf"

echo
echo "==> 2. Convert VCF to CSV (common + complex SV tables)"
python scripts/python/vcf_to_csv.py \
    "${OUTDIR}/${SAMPLE}.filtered.vcf" \
    --sample "${SAMPLE}" \
    --outdir "${OUTDIR}"

echo
echo "==> 3. Plot histograms (no_filter, >1kb, >10kb)"
python scripts/python/plot_histograms.py \
    "${OUTDIR}/${SAMPLE}_common_svs.csv" \
    --sample "${SAMPLE}" \
    --outdir "${HISTDIR}"

echo
echo "==> 4. Plot histograms restricted to chr9, chr17, chr18"
python scripts/python/plot_histograms.py \
    "${OUTDIR}/${SAMPLE}_common_svs.csv" \
    --sample "${SAMPLE}" \
    --chroms chr9,chr17,chr18 \
    --outdir "${HISTDIR}/filtered_chrs"

echo
echo "==> Done. See ${OUTDIR}/ for filtered VCF, CSVs, and histogram PNGs."
ls -la "${OUTDIR}"
echo
ls -la "${HISTDIR}"
