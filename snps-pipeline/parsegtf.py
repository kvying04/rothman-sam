import pandas as pd
from pathlib import Path
import json
import os

def gtfdf(gtfpath: Path) -> pd.DataFrame:
	gtfdf = pd.read_csv(gtfpath, sep='\t', comment='#')
	gtfdf.columns = ["chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
	
	return gtfdf


def getloc(gtf: pd.DataFrame) -> pd.DataFrame:
	locdf = pd.DataFrame()
	gtf = gtf[gtf['feature'] == 'gene']
	locdf['gene'] = gtf['attribute'].str.extract(r'gene_id "([^"]+)"', expand=False)
	locdf['chromosome'] = gtf['chromosome'].astype(str)
	locdf['start'] = gtf['start']
	locdf['end'] = gtf['end']

	return locdf

def convertchromosomes(locmappath: Path, chromosomemap: dict) -> None:
	locmap = pd.read_csv(locmappath, sep='\t')
	locmap['chromosome'] = locmap['chromosome'].str.split('.').str[0].map(chromosomemap).fillna(locmap['chromosome'])

	locmap.to_csv(locmappath, sep='\t', index=False)



if __name__ == "__main__":
	DIR = Path("/home/hanwenying/rothman")
	REFDIR = DIR / "ref"
	OUTDIR = DIR / "snps-pipeline/out"
	os.makedirs(OUTDIR, exist_ok=True)

	gtf = gtfdf(REFDIR/"ce11.gtf")
	locmap = getloc(gtf)
	locmap.to_csv(OUTDIR/"locmap.tsv", sep='\t', index=False)

	chromosomemap = pd.read_csv(REFDIR/"chromconversion.tsv", sep='\t')
	chromosomemap = chromosomemap.set_index('RefSeq')['Chrom'].to_dict()
	convertchromosomes(OUTDIR/"locmap.tsv", chromosomemap)