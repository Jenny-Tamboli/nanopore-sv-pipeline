# Nanopore Structural Variant Calling Pipeline

A reproducible pipeline for calling, filtering, classifying, and visualising structural variants (SVs) from Oxford Nanopore long-read sequencing data. Built around **NGMLR + Sniffles2**, with downstream classification into simple and complex SVs, chromosome-size-normalized SV-burden plots, and optional VEP gene-symbol annotation.

> **Note on data**: This repository contains code only. No proprietary or unpublished sequencing data is included. A synthetic VCF is provided under `tests/data/` so the downstream Python pipeline can be run end-to-end without access to real samples.

---

## Pipeline overview

```
FASTQ ──► NGMLR ──► BAM ──► sort/index ──► Sniffles2 ──► VCF
                                                          │
                                                          ▼
                                              filter (PASS + PRECISE,
                                                      numeric autosomes)
                                                          │
                                                          ▼
                                              classify + 1 kb filter
                                            ┌─────────────┴─────────────┐
                                            ▼                           ▼
                                    simple_svs.csv             complex_svs.csv
                                  (DEL/INS/DUP/INV,                  (BND)
                                   |SIZE| ≥ 1 kb)                      │
                                            │                           │
                          ┌─────────────────┼────────────┐              │
                          ▼                 ▼            ▼              ▼
                   SV burden plot      VEP annotate    Circa      VEP annotate
                  (normalized by      → gene symbols  (Circos)   → gene symbols
                   chrom size)
```

| Stage | Tool | Script |
|-------|------|--------|
| Alignment | NGMLR | `scripts/slurm/01_align_ngmlr.sh` |
| Sort + index | samtools | `scripts/slurm/02_sort_index.sh` |
| SV calling | Sniffles2 | `scripts/slurm/03_sniffles.sh` |
| VCF filtering | text | `scripts/python/filter_vcf.py` |
| Classify + 1 kb filter | pandas | `scripts/python/vcf_to_csv.py` |
| Burden plot | matplotlib | `scripts/python/plot_histograms.py` |
| VEP annotation merge | pandas | `scripts/python/annotate_vep.py` |

---

## Filtering and classification

**Stage 1 — VCF filter** (`filter_vcf.py`):
- Keep only `FILTER == PASS`
- Keep only records with `PRECISE` in INFO (drop IMPRECISE breakpoints)
- Numeric autosomes only — drop chrX, chrY, chrM, decoy/random contigs

**Stage 2 — classify and size-filter** (`vcf_to_csv.py`):
- `<sample>_simple_svs.csv` → DEL, INS, DUP, INV with `|SIZE| ≥ 1 kb`
- `<sample>_complex_svs.csv` → BND (no size filter — BNDs have no intrinsic length)

The 1 kb threshold is applied only to simple SVs because Sniffles reports SVLEN/END for those types. BNDs describe two breakend coordinates, so the same threshold doesn't apply. The `--min-size` flag in `vcf_to_csv.py` overrides the default 1000 bp.

---

## Output CSVs

**`<sample>_simple_svs.csv`** (one row per simple SV ≥ 1 kb)
| CHROM | START | STOP | SIZE | SVTYPE | STRAND | SIZE_kb |

`SIZE` preserves the signed SVLEN convention from Sniffles (negative for DEL). `SIZE_kb = SIZE / 1000`, also signed.

**`<sample>_complex_svs.csv`** (one row per BND)
| CHROM1 | START1 | CHROM2 | START2 | SVTYPE |

Both breakend ends are required to be on numeric autosomes; rows where either end is on chrX/Y/M/decoys are dropped.

---

## SV burden plot

`plot_histograms.py` produces a stacked horizontal bar chart, one bar per autosome, segmented by SVTYPE and capped by an "Other (unaffected)" remainder. Each segment width is the fraction of the chromosome affected:

    affected_fraction[svtype, chrom] = sum(|SIZE| for that svtype on that chrom)
                                       / chromosome_size

Chromosome sizes are parsed directly from the input VCF's `##contig=<ID=...,length=...>` header lines — no hardcoded reference assumptions. This normalization makes burden comparable across chromosomes regardless of their physical size.

Colours match the Circos plots (DEL red, INS yellow, DUP burgundy, INV green) so the two visualisations read consistently.

---

## VEP annotation (optional)

`annotate_vep.py` merges gene symbols from a VEP-annotated VCF onto either the simple- or complex-SVs CSV.

Intended workflow:
1. Upload `<sample>.filtered.vcf` to the [Ensembl VEP web server](https://www.ensembl.org/vep) with default parameters
2. Download the annotated VCF
3. Run `annotate_vep.py` once per CSV (simple / complex):

```bash
python scripts/python/annotate_vep.py \
    --sv-csv  results/sample/sample_simple_svs.csv \
    --vep-vcf results/sample/sample_vep.vcf \
    --mode    simple \
    --out     results/sample/sample_simple_svs_vep.csv
```

Join keys:
- **simple mode**: `(CHROM, START, STOP)` → inner join
- **complex mode**: `(CHROM1, START1)` → inner join, plus a `MATE_AGREES` column confirming the BND mate from VEP matches `CHROM2, START2` in your CSV

Both modes restrict to numeric autosomes, normalise the `chr` prefix convention to match your CSV, and use the Ensembl web server's default CSQ field order (overridable with `--csq-keys`).

---

## Circos plots

The CSVs are designed to feed directly into [Circa](http://omgenomics.com/circa/) or a Circos config. See `docs/circos_setup.md` for the track design:
- **Rectangle track** — SV span within a chromosome (from `simple_svs.csv`)
- **Scatter track** — SVTYPE ratio per chromosome
- **Connection track** — BND links between chromosomes (from `complex_svs.csv`)
- Colour code: DEL red, INS yellow, DUP burgundy, INV green

---

## Repository layout

```
nanopore-sv-pipeline/
├── README.md
├── LICENSE
├── environment.yml
├── requirements.txt
├── .gitignore
├── scripts/
│   ├── slurm/
│   │   ├── 01_align_ngmlr.sh
│   │   ├── 02_sort_index.sh
│   │   ├── 03_sniffles.sh
│   │   └── submit_all.sh
│   └── python/
│       ├── filter_vcf.py
│       ├── vcf_to_csv.py
│       ├── plot_histograms.py
│       └── annotate_vep.py
├── tests/
│   ├── data/example.vcf
│   └── test_pipeline.sh
├── docs/
│   ├── pipeline_overview.md
│   └── circos_setup.md
└── results/                # gitignored
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

### 2. Smoke test on synthetic data

```bash
bash tests/test_pipeline.sh
```

Skips alignment + SV calling. Runs filter → CSV → burden plot on `tests/data/example.vcf`. Produces filtered VCF, two CSVs, and a burden PNG in `results/test_sample/`.

### 3. Full pipeline on real data (SLURM)

```bash
bash scripts/slurm/submit_all.sh sample_name path/to/reads.fastq.gz path/to/reference.fa
```

Or run the same scripts with `bash` (not `sbatch`) outside a SLURM environment — the `#SBATCH` headers are inert. See `docs/pipeline_overview.md`.

### 4. Downstream

```bash
# Filter
python scripts/python/filter_vcf.py results/sample/sample.sniffles.vcf -o results/sample/sample.filtered.vcf

# Classify + 1 kb filter
python scripts/python/vcf_to_csv.py results/sample/sample.filtered.vcf --sample sample --outdir results/sample

# Burden plot
python scripts/python/plot_histograms.py results/sample/sample_simple_svs.csv \
    --vcf results/sample/sample.filtered.vcf --sample sample --outdir results/sample

# Optional: annotate with VEP output
python scripts/python/annotate_vep.py --sv-csv results/sample/sample_simple_svs.csv \
    --vep-vcf results/sample/sample_vep.vcf --mode simple --out results/sample/sample_simple_svs_vep.csv
```

---

## Citation

If you use this pipeline, please cite the underlying tools:
- **NGMLR / Sniffles**: Sedlazeck et al., *Nature Methods* (2018)
- **Sniffles2**: Smolka et al., *Nature Biotechnology* (2024)
- **samtools**: Danecek et al., *GigaScience* (2021)
- **VEP**: McLaren et al., *Genome Biology* (2016)

---

## License

MIT — see `LICENSE`.
