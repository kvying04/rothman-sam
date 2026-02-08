import parsegtf
import subprocess
import pandas as pd
from pathlib import Path

def removeheader(vcfpath: Path, outpath: Path) -> None:
	with open(vcfpath, 'r') as vcf, open(outpath, 'w') as cleanvcf:
		for line in vcf:
			if not line.startswith("##"):
				cleanvcf.write(line)

	return
	

def filterisolates(cleanpath: Path, isolates: list, filteredpath: Path) -> None:
	if filteredpath == None:
		filteredpath = cleanpath
	
	id = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

	cols = id + isolates

	filtereddf = pd.read_csv(cleanpath, usecols=cols, sep='\t')
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
	removeheader(vcfpath, cleanpath)
	print("\tFiltering isolates...")
	filterisolates(cleanpath, isolates, filteredpath)

	vcf = pd.read_csv(filteredpath, sep='\t')

	return vcf


if __name__ == "__main__":
	DIR = Path("/home/hanwenying/rothman")
	REFDIR = DIR / "ref"
	VCFPATH = REFDIR / "2025hardfilter.vcf"
	OUTDIR = DIR / "snps-pipeline/out"
	OUTPATH = OUTDIR / "2025clean.vcf"

	keep_strains = ['N2', 'LKC34', 'ED3017', 'JU775', 'MY16', 'MY23', 'JT11398', 'CB4856']
	loss_strains = ['CX11314']
	isolates = keep_strains + loss_strains

	vcfdf = loadvcf(VCFPATH, OUTPATH, isolates, None)
