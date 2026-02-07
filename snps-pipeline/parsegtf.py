import pandas as pd
from pathlib import Path
import os

def gtfdf(gtfpath: Path) -> pd.DataFrame:
	gtfdf = pd.read_csv(gtfpath, sep='\t', comment='#')
	gtfdf.columns = ["chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
	
	return gtfdf


def getloc(gtf: pd.DataFrame) -> pd.DataFrame:
	locdf = pd.DataFrame()

	locdf['gene'] = gtf['attribute'].str.extract(r'gene_id "([^"]+)"', expand=False)
	# location = gtf['start'].str.cat(gtf['end'].astype(str), sep="..")
	locdf['loc'] = gtf['chromosome'].astype(str) + ":" + gtf['start'].astype(str) + ".." + gtf['end'].astype(str)

	return locdf


if __name__ == "__main__":
	DIR = Path("/home/hanwenying/rothman")
	REFDIR = DIR / "ref"
	OUTDIR = DIR / "snps-pipeline/out"
	os.makedirs(OUTDIR, exist_ok=True)

	gtf = gtfdf(REFDIR/"ce11.gtf")
	locmap = getloc(gtf)
	locmap.to_csv(OUTDIR/"locmap.tsv", sep='\t', index=False)