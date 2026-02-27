from pathlib import Path
import os
import sys
import pandas as pd


ROOT = Path("/home/hanwenying/rothman-sam")

# importpath = ROOT / "snps-pipeline"
# sys.path.append(importpath)

import getsnps
import snpstogene

W76DIR = ROOT / "w76"
OUTDIR = W76DIR / "out"
LOCMAPPATH = "/home/hanwenying/rothman-sam/snps-pipeline/out/locmap.tsv"

conversionmap = pd.read_csv('/home/hanwenying/rothman-sam/ref/chromconversion.tsv', sep='\t')
# conversionmap = pd.Series(conversionmap['RefSeq'].values, index=conversionmap['Chrom']).to_dict()



vcf = getsnps.loadvcf(OUTDIR / "merged.vcf", OUTDIR / "merged_mod.vcf", ['A1', 'B2', 'C3', 'D4', 'E5', 'F6', 'G7', 'H8', 'I9', 'J10', 'K11', 'L12', 'M13', 'N14', 'O15', 'P16', 'Q17'], OUTDIR / "merged_mod.vcf")
isnps = getsnps.snpsofinterest(vcf, ['P16', 'A1', 'B2', 'C3', 'D4', 'E5', 'F6', 'G7', 'H8'], ['Q17', 'I9', 'J10', 'K11', 'L12', 'M13', 'N14', 'O15'])

c_dict = dict(zip(conversionmap['RefSeq'].str.strip(), conversionmap['Chrom']))

# 2. Process the isnps column
isnps['#CHROM'] = (
    isnps['#CHROM']
    .str.split('.')        # Split at the dot
    .str[0]                # Take the first part (NC_003279)
    .str.strip()           # Remove any hidden spaces
    .map(c_dict)           # Apply the map
    .fillna(isnps['#CHROM']) # Keep original if no match found
)

isnps.to_csv(OUTDIR/"interesting_w76.tsv", sep='\t', index=False)

keep = ['P16', 'A1', 'B2', 'C3', 'D4', 'E5', 'F6', 'G7', 'H8']
loss = ['Q17', 'I9', 'J10', 'K11', 'L12', 'M13', 'N14', 'O15']

snpgenes = snpstogene.snpstogene(OUTDIR / "interesting_w76.tsv", LOCMAPPATH, keep + loss)
snpgenes.to_csv(OUTDIR / "snpgene_w76.tsv", sep='\t', index=False)