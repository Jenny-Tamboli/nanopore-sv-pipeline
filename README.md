# Nanopore Structural Variant Calling Pipeline

A reproducible pipeline for calling, filtering, and visualising structural variants (SVs) from Oxford Nanopore long-read sequencing data. Built around **NGMLR + Sniffles2** with downstream filtering, VCF→CSV conversion, and histogram/Circos-ready outputs.

> **Note on data**: This repository contains code only. No proprietary or unpublished sequencing data is included. A small synthetic VCF is provided under `tests/data/` so all downstream scripts can be run end-to-end without access to real samples. For a realistic demo, see [Test data](#test-data) below for publicly available Nanopore datasets (e.g. GIAB HG002).

---

## Pipeline overview

```
FASTQ ──► NGMLR ──► BAM ──► sort/index ──► Sniffles2 ──► VCF
                                                          │
                                                          ▼
                                              filter (PASS + PRECISE)
                                                          │
                                                          ▼
                                              VCF → CSV (simple + complex SVs)
                                                          │
                                            ┌─────────────┴─────────────┐
                                            ▼                           ▼
                                    Histograms (matplotlib)     Circos input CSV
```

| Stage | Tool | Script |
|-------|------|--------|
| Alignment | NGMLR | `scripts/slurm/01_align_ngmlr.sh` |
| Sort + index | samtools | `scripts/slurm/02_sort_index.sh` |
| SV calling | Sniffles2 | `scripts/slurm/03_sniffles.sh` |
| VCF filtering | bcftools | `scripts/python/filter_vcf.py` |
| VCF → CSV | pysam / custom | `scripts/python/vcf_to_csv.py` |
| Histograms | matplotlib | `scripts/python/plot_histograms.py` |

---

## Repository layout

```
nanopore-sv-pipeline/
├── README.md
├── LICENSE
├── environment.yml          # conda environment
├── requirements.txt         # pip alternative
├── .gitignore
├── scripts/
│   ├── slurm/               # HPC job submission scripts
│   │   ├── 01_align_ngmlr.sh
│   │   ├── 02_sort_index.sh
│   │   ├── 03_sniffles.sh
│   │   └── submit_all.sh
│   └── python/              # Analysis scripts
│       ├── filter_vcf.py
│       ├── vcf_to_csv.py
│       └── plot_histograms.py
├── tests/
│   ├── data/                # Small synthetic VCFs for testing
│   │   └── example.vcf
│   └── test_pipeline.sh     # End-to-end demo on test data
├── docs/
│   ├── pipeline_overview.md
│   └── circos_setup.md
└── results/                 # Output directory (gitignored)
```

---

## Quick start

### 1. Install dependencies

```bash
conda env create -f environment.yml
conda activate nanopore-sv
```

Or with pip (Python parts only; you'll still need NGMLR, samtools, Sniffles2 installed separately):

```bash
pip install -r requirements.txt
```

### 2. Run the demo on synthetic data

```bash
bash tests/test_pipeline.sh
```

This skips alignment/SV calling and runs the Python pipeline on a synthetic VCF, producing filtered CSVs and histograms in `results/`.

### 3. Run the full pipeline on real data (SLURM)

Edit paths in `scripts/slurm/submit_all.sh` and submit:

```bash
bash scripts/slurm/submit_all.sh sample_name /path/to/reads.fastq.gz /path/to/reference.fa
```

---

## Filtering criteria

Following the report's spec, VCFs are filtered to retain only:
- Rows with `FILTER == PASS`
- Rows where `INFO` contains `PRECISE` (precise breakpoints; imprecise calls excluded)

Standard chromosomes only (`chr1`–`chr22`, `chrX`, `chrY`, `chrM`). Non-standard contigs like `chr18_gl000207_random` are dropped.

---

## Output CSVs

Two CSVs per sample:

**`<sample>_simple_svs.csv`** — DEL, INS, DUP, INV
| chrom | start | stop | size | svtype | strand |

**`<sample>_complex_svs.csv`** — BND (translocations / complex rearrangements)
| chrom1 | start1 | end1 | chrom2 | start2 | end2 | svtype |

> BNDs do not have an inherent length — they describe two breakend coordinates on (potentially) different chromosomes. Size-based filtering (`>1kb`, `>10kb`) is therefore applied only to simple SVs, not BNDs.

---

## Histograms

`plot_histograms.py` generates SV-size distributions per sample at three cutoffs:
- `no_filter`
- `size > 1 kb`
- `size > 10 kb`

One PNG per (sample, cutoff). Optional `--chroms chr9,chr17,chr18` restricts to a subset.

---

## Circos plots

CSVs are designed to feed directly into [Circa](http://omgenomics.com/circa/) or a Circos config. See `docs/circos_setup.md` for the track design:
- **Rectangle track** — SV span within a chromosome
- **Scatter track** — SVTYPE ratio per chromosome
- **Connection track** — BND links between chromosomes
- Colour code: DEL red, INS yellow, DUP burgundy, INV green

Use `no_filter` CSVs for Circos (so BNDs are included).

---

## Test data

For a realistic run without your own samples, download a small public Nanopore dataset:

- **nf-core/test-datasets (nanoseq branch)** — small ONT test BAM/FASTQ, fast to run.
- **GIAB HG002 ONT** — full reference sample, larger downloads. See `https://github.com/genome-in-a-bottle/giab_data_indexes`.
- **Sniffles2 GitHub** — ships a tiny demo dataset in `test/`.

---

## Citation

If you use this pipeline, please cite the underlying tools:
- **NGMLR / Sniffles**: Sedlazeck et al., *Nature Methods* (2018)
- **Sniffles2**: Smolka et al., *Nature Biotechnology* (2024)
- **samtools**: Danecek et al., *GigaScience* (2021)

---

## License

MIT — see `LICENSE`.
