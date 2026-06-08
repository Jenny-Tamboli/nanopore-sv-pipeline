import pandas as pd
import os
import numpy as np

sniffles = pd.read_csv('/Users/jennytamboli/Documents/carcinoma_normal_sniffles_filter/simple_csv/1kb_filter/hpde_1kb.csv')
vep = pd.read_csv('/Users/jennytamboli/Documents/carcinoma_normal_sniffles_filter/VEP/CSV/HPDE_vep_VEP.csv')
vep.rename(columns={'POS':'START', 'INFO_END': 'STOP'}, inplace=True)
vep['CHROM'] = vep['CHROM'].str.replace('chr', '')

# Filter out non-numeric chromosomes
vep = vep[vep['CHROM'].apply(lambda x: x.isnumeric())]

# Convert chromosome labels to integers (example: "1" to 1)
vep.loc[:, 'CHROM'] = vep['CHROM'].astype(int)

vep_ext = vep.loc[:,['CHROM', 'START', 'STOP', 'CSQ_SYMBOL']]

# Only include rows with matching keys in both DataFrames
annotate_commonSVs = pd.merge(sniffles, vep_ext, how='inner', left_on=['START', 'STOP'], right_on=['START', 'STOP']) 

df = annotate_commonSVs.drop(columns=['Unnamed: 0', 'CHROM_y']).rename(columns={'CHROM_x': 'CHROM'})

print(df)

df.to_csv('/Users/jennytamboli/Documents/carcinoma_normal_sniffles_filter/VEP/simple_variants_1kb/HPDE_1kb_vep.csv')

