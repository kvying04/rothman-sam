import pandas as pd
from pathlib import Path

if __name__ == "__main__":
	DIR = Path("/home/hanwenying/rothman-sam")
	OUTDIR = DIR / "snps-pipeline/out"
	# SNPGENEPATH = OUTDIR / "snpgene_w76.tsv"
	SNPGENEPATH = '/home/hanwenying/rothman-sam/w76/out/snpgene_w76.tsv'

	# genes = set('''C17H12.8
	# 		   C44B11.6
	# 		   D1086.7
	# 		   F07G6.10
	# 		   F14D2.13a
	# 		   F14D2.4a
	# 		   F30B5.4b
	# 		   F43E2.5
	# 		   F54D10.2
	# 		   F59A1.8
	# 		   M03C11.2
	# 		   M03C11.2
	# 		   R06B10.3
	# 		   R06B9.3
	# 		   R06B9.3
	# 		   R07B7.2b
	# 		   T05A10.3
	# 		   T05G11.1a
	# 		   T23B12.14
	# 		   W02D9.5
	# 		   Y17D7B.4
	# 		   Y43F8C.2
	# 		   Y65B4A.11
	# 		   ZK6.7c'''.split('\n'))
	with open(DIR / 'snps-pipeline/genelist.txt', 'r') as f:
		genes = [line.strip() for line in f]
	genes_fixed = []
	for gene in genes:
		if len(gene.split('.')) > 2:
			gene = '.'.join(gene.split('.')[:2])
		genes_fixed.append(gene)
	
	# genes = set(g.strip() for g in genes)
	
	snpgenes = pd.read_csv(SNPGENEPATH, sep='\t')
	snpgenes['Gene'] = snpgenes['Gene'].str.strip()
	intersections = snpgenes[snpgenes['Gene'].str.split('_').str[1].isin(genes_fixed)]

	intersections.to_csv(OUTDIR / "intersections_w76_2.tsv", sep='\t', index=False)