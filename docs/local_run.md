# Running locally (no SLURM)

If you don't have access to a SLURM cluster, the SLURM `#SBATCH` headers in
the scripts under `scripts/slurm/` are simply ignored when the scripts are
executed with `bash`. They work as ordinary shell scripts.

## Tier `tiny` — runs in minutes on a laptop

```bash
# 1. Get the data (small, public, NA12878 EDIL3 region)
bash scripts/download_test_data.sh tiny

# 2. Run alignment + sort/index + SV calling locally
bash scripts/slurm/01_align_ngmlr.sh tiny_demo \
     data/tiny/reads.fastq.gz \
     data/tiny/reference.fa

bash scripts/slurm/02_sort_index.sh tiny_demo

bash scripts/slurm/03_sniffles.sh tiny_demo data/tiny/reference.fa

# 3. Filter, convert, plot
python scripts/python/filter_vcf.py \
       results/tiny_demo/tiny_demo.sniffles.vcf \
       -o results/tiny_demo/tiny_demo.filtered.vcf

python scripts/python/vcf_to_csv.py \
       results/tiny_demo/tiny_demo.filtered.vcf \
       --sample tiny_demo --outdir results/tiny_demo

python scripts/python/plot_histograms.py \
       results/tiny_demo/tiny_demo_common_svs.csv \
       --sample tiny_demo --outdir results/tiny_demo/histplots
```

Because EDIL3 is a single gene, expect only a handful of SV calls — sometimes
none. The point of this tier is to prove every step of the plumbing works.

## Tier `hg002` — realistic chr22 demo

This tier uses pre-aligned reads from GIAB, so you can **skip steps 1 and 2**
and go straight to Sniffles2:

```bash
# 1. Get chr22 of HG002 (streamed, no 175 GB download)
bash scripts/download_test_data.sh hg002

# 2. Stage the pre-aligned BAM as if our pipeline produced it
mkdir -p results/hg002_chr22
cp data/hg002/HG002_chr22.bam     results/hg002_chr22/hg002_chr22.sorted.bam
cp data/hg002/HG002_chr22.bam.bai results/hg002_chr22/hg002_chr22.sorted.bam.bai

# 3. Call SVs
bash scripts/slurm/03_sniffles.sh hg002_chr22 data/hg002/reference.fa

# 4. Filter, convert, plot — same as above with the new sample name
python scripts/python/filter_vcf.py \
       results/hg002_chr22/hg002_chr22.sniffles.vcf \
       -o results/hg002_chr22/hg002_chr22.filtered.vcf

python scripts/python/vcf_to_csv.py \
       results/hg002_chr22/hg002_chr22.filtered.vcf \
       --sample hg002_chr22 --outdir results/hg002_chr22

python scripts/python/plot_histograms.py \
       results/hg002_chr22/hg002_chr22_common_svs.csv \
       --sample hg002_chr22 \
       --outdir results/hg002_chr22/histplots
```

Then load `hg002_chr22_common_svs.csv` and `hg002_chr22_complex_svs.csv` into
Circa (see `circos_setup.md`) for the Circos view.

## Running just the synthetic demo (no downloads)

For a zero-dependency smoke test that exercises only the Python parts:

```bash
bash tests/test_pipeline.sh
```

This skips alignment/SV calling and runs the filter → CSV → histogram chain
on the synthetic VCF in `tests/data/`.
