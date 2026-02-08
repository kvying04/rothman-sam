import pandas as pd
import pyranges as pr
from pathlib import Path

def snpstogene(snpspath: Path, locmappath: Path, isolates: list) -> pd.DataFrame:
	snps = pd.read_csv(snpspath, sep='\t')
	locmap = pd.read_csv(locmappath, sep='\t')

	snps = snps.loc[:, ['#CHROM', 'POS'] + isolates + [f"{isolate}_K" for isolate in isolates]]
	snps.columns = ['Chromosome', 'Start'] + isolates + [f"{isolate}_K" for isolate in isolates]

	locmap = locmap.loc[:, ['gene', 'chromosome', 'start', 'end']]
	locmap.columns = ['Gene', 'Chromosome', 'Start', 'End']

	snpspr = pr.PyRanges(snps.assign(End=snps['Start'] + 1))
	locmappr = pr.PyRanges(locmap)

	res = snpspr.join(locmappr)
	snpgenes = res.df
	snpgenes.drop(columns=['End'], inplace=True)
	snpgenes.columns = ['Chr', 'SNPPos'] + isolates + [f"{isolate}_K" for isolate in isolates] + ['Gene', 'Start', 'End']
	snpgenes = snpgenes[['Chr', 'SNPPos', 'Gene', 'Start', 'End'] + [f"{isolate}_K" for isolate in isolates] + isolates]
	snpgenes.drop_duplicates()

	return snpgenes

if __name__ == "__main__":
	DIR = Path("/home/hanwenying/rothman")
	OUTDIR = DIR / "snps-pipeline/out"
	SNPSPATH = OUTDIR / "interestingsnps.tsv"
	LOCMAPPATH = OUTDIR / "locmap.tsv"

	keep_strains = ['N2', 'LKC34', 'ED3017', 'JU775', 'MY16', 'MY23', 'JT11398', 'CB4856']
	loss_strains = ['CX11314']
	isolates = keep_strains + loss_strains
	
	snpgenes = snpstogene(SNPSPATH, LOCMAPPATH, isolates)
	snpgenes.to_csv(OUTDIR/"snpgene.tsv", sep='\t', index=False)