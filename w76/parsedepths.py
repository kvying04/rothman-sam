import pandas as pd
from pathlib import Path

ROOT = Path('/home/hanwenying')
DEPTHS = ROOT / "rothman-sam/w76/out/depths"

genes = ['all']
bams = ['A1', 'B2', 'C3', 'D4', 'E5', 'F6', 'G7', 'H8', 'I9', 'J10', 'K11', 'L12', 'M13', 'N14', 'O15', 'P16', 'Q17']


for gene in genes:
	for bam in bams:
		path = DEPTHS / f"{bam}-{gene}.tsv"
		df = pd.read_csv(path, sep='\t', header=None)
		df.columns = ['CHROM', 'POS', 'DEPTH']
		notread = df[df['DEPTH'] == 0]
		print(f"{bam}-{gene}: {len(notread)}/{len(df)}")