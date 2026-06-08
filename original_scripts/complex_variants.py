import pandas as pd
import os
import vcfpy
import numpy as np

sniffles = pd.read_csv('/Users/jennytamboli/Documents/carcinoma_normal_sniffles_filter/complex_csv/Panc1.csv')
vep = pd.read_csv('/Users/jennytamboli/Documents/carcinoma_normal_sniffles_filter/VEP/Panc1_vep_VEP.csv')
vep['CHROM'] = vep['CHROM'].str.replace('chr', '')

# Filter out non-numeric chromosomes (chrX, chrY, chrM)
vep = vep[vep['CHROM'].apply(lambda x: x.isnumeric())]
vep.loc[:, 'CHROM'] = vep['CHROM'].astype(int)

# Extract specific chromosomes and BND type SVs
selective_chr = vep[(vep['CHROM'].isin([9, 15, 17, 18])) & (vep['INFO_SVTYPE'] == 'BND')]
circ = selective_chr.loc[:, ['CHROM','POS','CSQ_Allele','CSQ_SYMBOL','INFO_CHR2']].copy()
circ['CSQ_Allele'] = circ['CSQ_Allele'].replace('.N', np.nan)


# Function to extract the numeric chromosomal location
def extract_numeric_location(allele):
    if pd.isna(allele):
        return None
    try:
        # Split first by '[' and take the second part, then split by ':' and take the second part (the number)
        return allele.split('[')[1].split(':')[1].split('[')[0]
    except IndexError:
        return None

# Apply the function to the CSQ_Allele column
circ['POS2'] = circ['CSQ_Allele'].apply(extract_numeric_location)
circ.drop(columns=['CSQ_Allele'], inplace=True)

circ['INFO_CHR2'] = circ['INFO_CHR2'].str.replace('chr', '')
circ = circ[circ['INFO_CHR2'].apply(lambda x: x.isnumeric())]
circ['INFO_CHR2'] = circ['INFO_CHR2'].astype(int)

# manipulate converted csv file from sniffles with just BNDs

sniffles['CHROM1'] = sniffles['CHROM1'].str.replace('chr', '')
sniffles['CHROM2'] = sniffles['CHROM2'].str.replace('chr', '')
sniffles = sniffles[sniffles['CHROM1'].apply(lambda x: x.isnumeric())]
sniffles = sniffles[sniffles['CHROM2'].apply(lambda x: x.isnumeric())]
sniffles['CHROM1'] = sniffles['CHROM1'].astype(int)
sniffles['CHROM2'] = sniffles['CHROM2'].astype(int)

merge_snif_vep = pd.merge(sniffles, circ, how='inner', left_on=['CHROM1', 'START1'], right_on=['CHROM', 'POS'])
merge_snif_vep.drop(columns=['CHROM', 'POS', 'INFO_CHR2', 'POS2'], inplace=True)
print(merge_snif_vep)
merge_snif_vep.to_csv('/Users/jennytamboli/Documents/carcinoma_normal_sniffles_filter/VEP/Panc1_bnd_vep.csv')
