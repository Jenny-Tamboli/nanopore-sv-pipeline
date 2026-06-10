#!/bin/bash
# -----------------------------------------------------------------------------
# Download publicly available Nanopore data for testing the pipeline.
#
# Three tiers, pick what fits your machine / use case:
#
#   tiny    nf-core/test-datasets nanoseq branch — NA12878 ONT reads aligned
#           to the EDIL3 gene region. A few MB total. Runs end-to-end in
#           minutes on a laptop. Best for CI and quick demos. Single gene, so
#           don't expect many SV calls but the pipeline mechanics all run.
#
#   hg002   GIAB HG002 ultralong ONT — stream just chr22 from the hosted
#           GRCh38 BAM via samtools (no full 175 GB download). Realistic SV
#           landscape across one chromosome. Requires samtools + ~3-5 GB
#           disk + ~30-60 min over a decent connection.
#
#   norris  Norris et al. 2016 — published Nanopore SV-detection dataset
#           on pancreatic cancer cell-line SVs (CDKN2A, SMAD4, etc.).
#           SRA SRP069199. Amplicon data with ~10 known SVs validated by
#           the original authors. Requires sra-toolkit + ngmlr + samtools.
#           Runs the full alignment+sort+SV-call pipeline (~1-2 hours).
#           Best for demonstrating the pipeline on real published data.
#
# Usage:
#   bash scripts/download_test_data.sh tiny
#   bash scripts/download_test_data.sh hg002
#   bash scripts/download_test_data.sh norris
#
# Output:
#   data/<tier>/        downloaded reads and reference
#   data/<tier>/README  provenance notes (sources, dates, citation)
# -----------------------------------------------------------------------------

set -euo pipefail

TIER="${1:-tiny}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"

# --- helpers ----------------------------------------------------------------
have() { command -v "$1" >/dev/null 2>&1; }

require() {
    local cmd="$1"
    if ! have "${cmd}"; then
        echo "ERROR: '${cmd}' is required but not found in PATH." >&2
        echo "  Install via conda (see environment.yml) or your package manager." >&2
        exit 1
    fi
}

download() {
    local url="$1" dest="$2"
    if [[ -s "${dest}" ]]; then
        echo "  -> already present, skipping: ${dest}"
        return 0
    fi
    echo "  -> ${url}"
    if have curl; then
        curl -fL --retry 3 --retry-delay 5 -o "${dest}" "${url}"
    elif have wget; then
        wget -c -O "${dest}" "${url}"
    else
        echo "ERROR: need curl or wget." >&2
        exit 1
    fi
}

# --- tier: tiny -------------------------------------------------------------
download_tiny() {
    # nf-core/test-datasets nanoseq branch — NA12878 ONT reads overlapping
    # the EDIL3 gene, plus a sub-setted GRCh38 EDIL3 reference. Hosted on
    # GitHub, so URLs are very stable.
    local OUTDIR="data/tiny"
    local BASE="https://raw.githubusercontent.com/nf-core/test-datasets/nanoseq"

    mkdir -p "${OUTDIR}"

    echo "[tiny] downloading EDIL3 reference (GRCh38)"
    download "${BASE}/reference/GRCh38_EDIL3.fa" "${OUTDIR}/reference.fa"

    echo "[tiny] downloading NA12878 ONT FASTQ (EDIL3 region)"
    download "${BASE}/fastq/nondemultiplexed/sample_nobc.fastq.gz" \
             "${OUTDIR}/reads.fastq.gz"

    cat > "${OUTDIR}/README.md" <<'EOF'
# Tier: tiny — nf-core nanoseq test data

NA12878 Nanopore reads overlapping the EDIL3 gene, plus a region-specific
GRCh38 reference. Sourced from the nf-core/test-datasets repo, nanoseq
branch. Public, no embargo.

Files:
  reference.fa       GRCh38 EDIL3 region
  reads.fastq.gz     NA12878 ONT reads (nondemultiplexed subset)

Citation / source:
  https://github.com/nf-core/test-datasets/tree/nanoseq
EOF

    echo
    echo "[tiny] done. Try:"
    echo "  bash scripts/slurm/submit_all.sh tiny_demo \\"
    echo "    ${OUTDIR}/reads.fastq.gz \\"
    echo "    ${OUTDIR}/reference.fa"
    echo
    echo "Or for a non-SLURM local run, see docs/local_run.md."
}

# --- tier: hg002 ------------------------------------------------------------
download_hg002() {
    # Stream just chr22 from the GIAB HG002 GRCh38 ultralong ONT BAM. This
    # avoids the full 175 GB download — samtools uses HTTP range requests
    # against the public NIST FTP mirror.
    require samtools
    require curl

    local OUTDIR="data/hg002"
    mkdir -p "${OUTDIR}"

    local BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/HG002_GRCh38_ONT-UL_UCSC_20200508.phased.bam"
    local CHR="chr22"

    echo "[hg002] downloading GRCh38 reference (no_alt analysis set)"
    download \
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz" \
        "${OUTDIR}/reference.fa.gz"

    if [[ ! -s "${OUTDIR}/reference.fa" ]]; then
        echo "[hg002] decompressing reference"
        gunzip -k "${OUTDIR}/reference.fa.gz"
    fi

    if [[ ! -s "${OUTDIR}/reference.fa.fai" ]]; then
        echo "[hg002] indexing reference"
        samtools faidx "${OUTDIR}/reference.fa"
    fi

    local OUT_BAM="${OUTDIR}/HG002_chr22.bam"
    if [[ ! -s "${OUT_BAM}" ]]; then
        echo "[hg002] streaming ${CHR} from GIAB HG002 BAM (this can take a while)"
        # samtools view supports remote URLs and a region; it pulls only the
        # byte ranges needed for that chromosome.
        samtools view -b -@ 4 "${BAM_URL}" "${CHR}" -o "${OUT_BAM}"
        samtools index -@ 4 "${OUT_BAM}"
    else
        echo "[hg002] ${OUT_BAM} already exists, skipping stream"
    fi

    cat > "${OUTDIR}/README.md" <<EOF
# Tier: hg002 — GIAB HG002 chr22 (ONT ultralong)

Reads from the GIAB Ashkenazi son HG002 (GM24385) ultralong Nanopore
release, restricted to ${CHR}. Aligned to GRCh38 by UCSC.

Files:
  reference.fa            GRCh38 no_alt analysis set (full)
  reference.fa.fai        samtools index
  HG002_chr22.bam         streamed chr22 alignments
  HG002_chr22.bam.bai     BAM index

These alignments already exist (i.e. step 1 — NGMLR — is not needed).
To call SVs:

  sbatch scripts/slurm/03_sniffles.sh hg002_chr22 \\
         ${OUTDIR}/reference.fa
  # (after copying/symlinking HG002_chr22.bam -> results/hg002_chr22/hg002_chr22.sorted.bam)

Citation / source:
  https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/
  README: README_ONT-UL_UCSC_HG002.md (in the same directory)
EOF

    echo
    echo "[hg002] done. The BAM is already sorted+indexed; you can skip"
    echo "        alignment + sort and go straight to Sniffles2:"
    echo
    echo "  mkdir -p results/hg002_chr22"
    echo "  cp ${OUT_BAM}     results/hg002_chr22/hg002_chr22.sorted.bam"
    echo "  cp ${OUT_BAM}.bai results/hg002_chr22/hg002_chr22.sorted.bam.bai"
    echo "  sbatch scripts/slurm/03_sniffles.sh hg002_chr22 ${OUTDIR}/reference.fa"
}

# --- tier: norris (Norris et al. 2016, SRP069199) ---------------------------
download_norris() {
    # Norris AL, Workman RE, Fan Y, Eshleman JR, Timp W. (2016)
    # "Nanopore sequencing detects structural variants in cancer."
    # Cancer Biol Ther. 17(3):246-253. DOI: 10.1080/15384047.2016.1139236
    # SRA: SRP069199 — ONT MinION 2D reads of pancreatic-cancer SV amplicons.
    require samtools
    require ngmlr
    require sniffles
    if ! have prefetch || ! have fasterq-dump; then
        echo "ERROR: sra-toolkit (prefetch, fasterq-dump) required." >&2
        echo "  Install with: conda install -c bioconda sra-tools" >&2
        exit 1
    fi

    local OUTDIR="data/norris2016"
    mkdir -p "${OUTDIR}"

    # Specific runs from SRP069199. SRR3061551 is a small representative run;
    # change SRA_RUNS to add more samples.
    local SRA_RUNS=("SRR3061551")

    # hg19 reference (paper used hg19 + BWA-MEM; we use NGMLR for parity with
    # the rest of the pipeline).
    local REF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
    local REF_GZ="${OUTDIR}/hg19.fa.gz"
    local REF_FA="${OUTDIR}/hg19.fa"

    echo "[norris] downloading hg19 reference (~900 MB compressed, ~3 GB decompressed)"
    if [[ ! -s "${REF_FA}" ]]; then
        download "${REF_URL}" "${REF_GZ}"
        echo "[norris] decompressing reference"
        gunzip -k "${REF_GZ}"
    else
        echo "  -> reference already present"
    fi

    if [[ ! -s "${REF_FA}.fai" ]]; then
        echo "[norris] indexing reference"
        samtools faidx "${REF_FA}"
    fi

    # Fetch each SRA run as FASTQ
    for SRR in "${SRA_RUNS[@]}"; do
        local FQ="${OUTDIR}/${SRR}.fastq.gz"
        if [[ -s "${FQ}" ]]; then
            echo "[norris] ${FQ} already present, skipping fetch"
            continue
        fi
        echo "[norris] fetching ${SRR} from SRA"
        ( cd "${OUTDIR}" && prefetch "${SRR}" && fasterq-dump --threads 4 "${SRR}" )
        gzip "${OUTDIR}/${SRR}.fastq" || true
        # Clean up the .sra cache
        rm -rf "${OUTDIR}/${SRR}/" 2>/dev/null || true
    done

    cat > "${OUTDIR}/README.md" <<EOF
# Tier: norris2016 — Norris et al. 2016 (SRP069199)

ONT MinION 2D reads of pancreatic-cancer structural variant amplicons,
covering known SVs in CDKN2A/p16, SMAD4/DPC4, and others.

Citation:
  Norris AL, Workman RE, Fan Y, Eshleman JR, Timp W.
  Nanopore sequencing detects structural variants in cancer.
  Cancer Biol Ther. 2016;17(3):246-253.
  DOI: 10.1080/15384047.2016.1139236
  SRA:  SRP069199

Files:
  hg19.fa                hg19 reference (matches the paper's alignment)
  hg19.fa.fai            samtools index
  SRR*.fastq.gz          ONT reads, one file per SRA run

Note: The paper used BWA-MEM with -x ont2d. We use NGMLR for consistency
with the rest of this pipeline, which is also Sniffles2-compatible. SV
calls will not be identical to those in the paper (Sniffles2 vs. custom
split-read analysis) but should recover the same large events.

To run the pipeline (local, no SLURM):

  SAMPLE=norris_demo
  bash scripts/slurm/01_align_ngmlr.sh \$SAMPLE \\
       ${OUTDIR}/SRR3061551.fastq.gz \\
       ${OUTDIR}/hg19.fa
  bash scripts/slurm/02_sort_index.sh \$SAMPLE
  bash scripts/slurm/03_sniffles.sh   \$SAMPLE ${OUTDIR}/hg19.fa

After the run, you'll have results/\$SAMPLE/\$SAMPLE.sniffles.vcf — a real
Sniffles2 VCF from published data. To use it as the repo's smoke-test VCF:

  cp results/\$SAMPLE/\$SAMPLE.sniffles.vcf tests/data/norris2016.vcf
  # Then edit tests/test_pipeline.sh to point at the new VCF.
EOF

    echo
    echo "[norris] download complete. Run the pipeline with:"
    echo "  bash scripts/slurm/01_align_ngmlr.sh norris_demo \\"
    echo "       ${OUTDIR}/SRR3061551.fastq.gz ${OUTDIR}/hg19.fa"
    echo "  bash scripts/slurm/02_sort_index.sh norris_demo"
    echo "  bash scripts/slurm/03_sniffles.sh   norris_demo ${OUTDIR}/hg19.fa"
    echo
    echo "See ${OUTDIR}/README.md for full instructions."
}

# --- dispatch ---------------------------------------------------------------
case "${TIER}" in
    tiny)   download_tiny ;;
    hg002)  download_hg002 ;;
    norris) download_norris ;;
    -h|--help|help)
        sed -n '2,33p' "$0"
        exit 0
        ;;
    *)
        echo "Unknown tier: '${TIER}'. Choose: tiny | hg002 | norris" >&2
        exit 2
        ;;
esac
