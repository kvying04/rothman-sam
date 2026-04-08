import pandas as pd
from pathlib import Path

def _removeheader(vcfpath: Path, outpath: Path) -> None:
	with open(vcfpath, 'r') as vcf, open(outpath, 'w') as cleanvcf:
		for line in vcf:
			if not line.startswith("##"):
				cleanvcf.write(line)

	return
	

def _filterisolates(cleanpath: Path, isolates: list, filteredpath: Path) -> None:
	if filteredpath == None:
		filteredpath = cleanpath
	
	id = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

	cols = id + isolates

	print("\t\tReading VCF into Pandas...")
	filtereddf = pd.read_csv(cleanpath, usecols=cols, sep='\t')
	print("\t\tWriting to CSV...")
	filtereddf.to_csv(filteredpath, sep='\t', index=False)
	# colindices = []
	
	# for col in cols:
	# 	try:
	# 		colindices.append(str(colnames.index(col) + 1))
	# 	except ValueError:
	# 		print(f"'{col}' is not in the VCF.")

	# colindices = ','.join(colindices)
	# print(colindices)

	# with open(filteredpath, "w") as f:
	# 	subprocess.run(["cut", "-d", "\t", "-f", colindices, cleanpath], stdout=f)

	return


def loadvcf(vcfpath: Path, cleanpath: Path, isolates: list, filteredpath: Path) -> pd.DataFrame:
	print("Remaking VCF...")
	print("\tRemoving header rows...")
	_removeheader(vcfpath, cleanpath)
	print("\tFiltering isolates...")
	_filterisolates(cleanpath, isolates, filteredpath)

	vcf = pd.read_csv(filteredpath, sep='\t')

	return vcf


def _snpsofinterest_col(vcfcol: pd.Series, phenotype: str) -> pd.Series:
	vcfcol = vcfcol.str.extract(r'(\d+)/(\d+)')
	if phenotype == 'loss':
		vcfcol_summable = vcfcol.fillna(1)
	elif phenotype == 'keep':
		vcfcol_summable = vcfcol.fillna(0)
	
	vcfcol = vcfcol.fillna(-1).astype(int).sum(axis=1)
	vcfcol = vcfcol.mask(vcfcol < 0, '.')

	vcfcol_summable = vcfcol_summable.astype(int).sum(axis=1)
	
	return vcfcol_summable, vcfcol


def snpsofinterest(vcf: pd.DataFrame, keep_strains: list, loss_strains: list) -> pd.DataFrame:
	interestingsnps = vcf.loc[:, '#CHROM':'FORMAT']
	snpfilter = pd.DataFrame()

	for isolate in keep_strains:
		interestingsnps[isolate] = vcf[isolate]
		snpfilter[isolate], interestingsnps[f"{isolate}_condition"] = _snpsofinterest_col(vcf[isolate], 'keep')

	for isolate in loss_strains:
		interestingsnps[isolate] = vcf[isolate]
		snpfilter[isolate], interestingsnps[f"{isolate}_condition"] = _snpsofinterest_col(vcf[isolate], 'loss')

	snpfilter['max_keep'] = snpfilter[keep_strains].max(axis=1)
	snpfilter['min_loss'] = snpfilter[loss_strains].min(axis=1)

	return interestingsnps[snpfilter['max_keep'] < snpfilter['min_loss']]


if __name__ == "__main__":
	DIR = Path("/home/hanwenying/rothman")
	REFDIR = DIR / "ref"
	VCFPATH = REFDIR / "2025hardfilter.vcf"
	OUTDIR = DIR / "snps-pipeline/out"
	OUTPATH = OUTDIR / "2025clean.vcf"

	keep_strains = ['N2', 'LKC34', 'ED3017', 'JU775', 'MY16', 'MY23', 'JT11398', 'CB4856', 'P16', 'A1', 'B2', 'C3', 'D4', 'E5', 'F6', 'G7', 'H8']
	loss_strains = ['CX11314', 'Q17', 'I9', 'J10', 'K11', 'L12', 'M13', 'N14', 'O15']
	isolates = keep_strains + loss_strains

	# vcfdf = loadvcf(VCFPATH, OUTPATH, isolates, None)
	vcfdf = pd.read_csv(OUTPATH, sep='\t')

	isnps = snpsofinterest(vcfdf, keep_strains, loss_strains)
	isnps.to_csv(OUTDIR/"interestingsnps.tsv", sep='\t', index=False)