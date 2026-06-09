#!/bin/bash
# -----------------------------------------------------------------------------
# End-to-end demo: filter -> CSV (1 kb) -> chromosome-normalized burden plot
# on the synthetic example VCF. Run from the repo root:
#   bash tests/test_pipeline.sh
# -----------------------------------------------------------------------------

set -euo pipefail

SAMPLE="test_sample"
OUTDIR="results/${SAMPLE}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"

mkdir -p "${OUTDIR}"

echo "==> 1. Filter VCF (PASS + PRECISE, numeric autosomes only)"
python scripts/python/filter_vcf.py \
    tests/data/example.vcf \
    -o "${OUTDIR}/${SAMPLE}.filtered.vcf"

echo
echo "==> 2. Convert VCF to CSVs (simple SVs with 1 kb filter + complex SVs)"
python scripts/python/vcf_to_csv.py \
    "${OUTDIR}/${SAMPLE}.filtered.vcf" \
    --sample "${SAMPLE}" \
    --outdir "${OUTDIR}"

echo
echo "==> 3. Plot per-chromosome SV burden (size-normalized)"
python scripts/python/plot_histograms.py \
    "${OUTDIR}/${SAMPLE}_simple_svs.csv" \
    --vcf "${OUTDIR}/${SAMPLE}.filtered.vcf" \
    --sample "${SAMPLE}" \
    --outdir "${OUTDIR}"

echo
echo "==> Done. Files in ${OUTDIR}/:"
ls -la "${OUTDIR}"
