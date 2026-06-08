import pandas as pd
import os
import vcfpy
import numpy as np


# Directory containing the VCF files
vcf_dir = '/Users/jennytamboli/Documents/carcinoma_normal_sniffles_filter/VEP/VCF'

# Directory where the converted CSV files will be saved
output_dir = '/Users/jennytamboli/Documents/carcinoma_normal_sniffles_filter/VEP/CSV'

# List all VCF files in the directory
#vcf_files = [f for f in os.listdir(vcf_dir) if f.endswith('.vcf')]
vcf_files = [f for f in os.listdir(vcf_dir) if f.startswith('HPDE')]

# Define the CSQ keys
csq_keys = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", 
            "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", 
            "Amino_acids", "Codons", "Existing_variation", "REF_ALLELE", "UPLOADED_ALLELE", "DISTANCE", 
            "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID", "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", 
            "APPRIS", "SIFT", "PolyPhen", "AF", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME",
            "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "TRANSCRIPTION_FACTORS"
           ]

# Iterate through each VCF file and process it
for vcf_file in vcf_files:
    # Construct full file path
    vcf_file_path = os.path.join(vcf_dir, vcf_file)
    
    # Create VCF reader
    reader = vcfpy.Reader(open(vcf_file_path, 'r'))
    
    data = []
    
    # Process each record in the VCF file
    for record in reader:
        chrom = record.CHROM
        pos = record.POS
        info = record.INFO

        row = {
            'CHROM': chrom,
            'POS': pos,
        }

        # Process CSQ field if it exists
        if 'CSQ' in info:
            csq_values = info['CSQ'][0].split('|')
            for csq_key, csq_value in zip(csq_keys, csq_values):
                row[f'CSQ_{csq_key}'] = csq_value

        # Add other INFO fields to the row dictionary
        for key, value in info.items():
            if key != 'CSQ':    # Skip CSQ since it's already processed
                row[f'INFO_{key}'] = value

        data.append(row)  # Add the row dictionary to the list

    # Convert the list of rows to a DataFrame
    df = pd.DataFrame(data)
    
    # Create output CSV file name
    output_csv_name = os.path.splitext(vcf_file)[0] + '_VEP.csv'
    output_csv_path = os.path.join(output_dir, output_csv_name)
    
    # Save the DataFrame to a CSV file
    df.to_csv(output_csv_path, index=False)
    
    print(f'Processed {vcf_file} and saved to {output_csv_path}')
