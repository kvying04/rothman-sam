import pandas as pd
from pathlib import Path
import numpy as np

ROOT = Path("/home/hanwenying")
INDELDIR = ROOT / "rothman-sam/w76/imsindel-out"
OUTDIR = ROOT / "rothman-sam/w76/out/indelouts"

lowgroup = ['A1', 'B2', 'C3']
highgroup = ['M13', 'N14', 'O15']

chroms = ['NC_003279.8', 'NC_003280.10', 'NC_003281.10', 'NC_003282.8', 'NC_003283.11', 'NC_003284.9', 'NC_001328.1']

# for each bam thing (e.g. A1): read in NC_003279.8.out; intersection for the lowgroup; set diff the union of highgroup
# telling "same indel": startbp within 50bp of each other

def compare_chrom_1(chrom: str, lowgroup: list, highgroup: list) -> pd.DataFrame: 
	# reading in indels
	lowindels = []
	highindels = []

	for bam in lowgroup:
		try:
			df = pd.read_csv(INDELDIR / f"{bam}.bam/{chrom}.out", sep='\t')
			df['source'] = bam
			lowindels.append(df)
		except FileNotFoundError:
			print(f"{bam}.bam/{chrom}.out does not exist; continuing...")
			pass
	for bam in highgroup:
		try:
			df = pd.read_csv(INDELDIR / f"{bam}.bam/{chrom}.out", sep='\t')
			df['source'] = bam
			highindels.append(df)
		except FileNotFoundError:
			print(f"{bam}.bam/{chrom}.out does not exist; continuing...")
			pass

	lowdf = pd.concat(lowindels, ignore_index=True)
	highdf = pd.concat(highindels, ignore_index=True)

	lowdf = lowdf[lowdf.groupby('sttpos')['sttpos'].transform('count') == (len(lowgroup))] # -5 has some
	lowdf['sttpos'] = lowdf['sttpos'].astype(np.int32)
	highstarts = list(set((highdf["sttpos"].astype(np.int32))))

	lowdf = lowdf[~lowdf["sttpos"].isin(highstarts)]

	return lowdf


for chrom in chroms:
	c = compare_chrom_1(chrom, lowgroup, highgroup)
	if len(c) > 0:
		print(chrom, len(c))
		c.to_csv(OUTDIR / f"{chrom}.tsv", sep='\t', index=False)