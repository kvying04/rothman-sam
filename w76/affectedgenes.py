import pandas as pd
from pathlib import Path

ROOT = Path("/home/hanwenying")
W76DIR = ROOT / "rothman-sam/w76"
REFDIR = ROOT / "rothman-sam/ref"

gtfdf = pd.read_csv(REFDIR / "ce11.gtf", sep='\t', skiprows=4, header=None)

gtfdf.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

gtfdf[gtfdf['feature'] == 'gene']
gtfdf['chrom'] = gtfdf['chrom'].str.split('.').str[0]

exindels = pd.read_csv(W76DIR / "imsindel-out/A1.bam/NC_003279.8.out", sep='\t')

exindels['chr'] = exindels['chr'].str.split('.').str[0]

def checkoverlap(genestart, geneend, indelstart, indelend):
	return max(genestart, indelstart) < min(geneend, indelend)

def searchindels(indeldf, gtfdf):
	cols = list(indeldf.columns)
	cols.append('affectedgene')
	outdf = pd.DataFrame(columns=cols)
	chr = exindels['chr'].unique()[0] # should only be len=1
	gtfdf = gtfdf[gtfdf['chrom'] == chr]
	for i, row in indeldf.iterrows():
		indel_start, indel_end = row['sttpos'], row['endpos']
		for i, generow in gtfdf.iterrows():
			gene_start, gene_end = generow['start'], generow['end']

			if checkoverlap(gene_start, gene_end, indel_start, indel_end):
				row['affectedgene'] = generow['attribute']
				outdf.loc[len(outdf)] = row
	
	return outdf


bams = ['A1.bam', 'B2.bam', 'C3.bam', 'D4.bam', 'E5.bam', 'F6.bam', 'G7.bam', 'H8.bam', 'I9.bam', 'J10.bam', 'K11.bam', 'L12.bam', 'M13.bam', 'N14.bam', 'O15.bam', 'P16.bam', 'Q17.bam']
chroms = ['NC_001328.1', 'NC_003279.8', 'NC_003280.10', 'NC_003281.10', 'NC_003282.8', 'NC_003283.11', 'NC_003284.9']

cols = list(exindels.columns)
cols.append('affectedgene')

def chromsrun(bamfolder):
	megaoutdf = pd.DataFrame(columns=cols)
	for chrom in chroms:
		print("\t"+chrom)
		indeldf = pd.read_csv(W76DIR / f"imsindel-out/{bamfolder}/{chrom}.out", sep='\t')
		indeldf['chr'] = indeldf['chr'].str.split('.').str[0]
		outdf = searchindels(indeldf, gtfdf)
		megaoutdf = pd.concat([megaoutdf, outdf], ignore_index=True)
	megaoutdf.to_csv(W76DIR / "affectedgenes" / f"{bamfolder}_affectedgenes.tsv", sep='\t', index=False)

# for bamfolder in bams:
# 	megaoutdf = pd.DataFrame(columns=cols)
# 	print(bamfolder)
# 	for chrom in chroms:
# 		print("\t"+chrom)
# 		indeldf = pd.read_csv(W76DIR / f"imsindel-out/{bamfolder}/{chrom}.out", sep='\t')
# 		indeldf['chr'] = indeldf['chr'].str.split('.').str[0]
# 		outdf = searchindels(indeldf, gtfdf)
# 		megaoutdf = pd.concat([megaoutdf, outdf], ignore_index=True)
# 	megaoutdf.to_csv(W76DIR / "affectedgenes" / f"{bamfolder}_affectedgenes.tsv", sep='\t', index=False)

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

if __name__ == "__main__":
	with ProcessPoolExecutor(max_workers=16) as exec:
		out = exec.map(chromsrun, bams)