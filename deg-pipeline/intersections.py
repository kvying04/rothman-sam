import pandas as pd
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import time


def _targetstrains(tpms: pd.DataFrame, loss_strains: list, keep_strains: list) -> pd.DataFrame:
	cols = ['transcript', 'N2'] + loss_strains + keep_strains
	tpms.columns = [col.split('_')[0] for col in tpms.columns]
	tpms = tpms[cols]
	# tpms.sort_values(by='transcript', inplace=True)
	tpms = tpms.set_index('transcript')
	tpms.columns = [col + f'_{(i%(sum(tpms.columns == col)))+1}' for i, col in enumerate(tpms.columns)]

	return tpms


def _metacount(tpms: pd.DataFrame, strain: str) -> tuple[pd.DataFrame, pd.DataFrame]:
	counts = tpms.loc[:, tpms.columns.str.startswith(('N2',strain))].round().astype(int)
	metadata = pd.DataFrame(index=counts.columns)
	metadata['condition'] = [exp.split('_')[0] for exp in metadata.index]

	return counts, metadata


def _dds_strain(counts: pd.DataFrame, metadata: pd.DataFrame, design_factors: str='condition', refit_cooks: bool=True) -> DeseqDataSet:
	dds = DeseqDataSet(	counts=counts.T,
						metadata=metadata,
						design_factors=design_factors,
						refit_cooks=refit_cooks,
						quiet=True,
						n_cpus=1)
	
	dds.deseq2()

	return dds


def _ds_strain(dds: DeseqDataSet, strain: str) -> pd.DataFrame: # take in DeseqStats params?
	ds = DeseqStats(dds, contrast=["condition", strain, "N2"], quiet=True, n_cpus=1)
	ds.summary()
	deg = ds.results_df.copy()
	
	return deg


def _rundeseq(strain: str, metacounts: tuple[pd.DataFrame, pd.DataFrame], design_factors: str='condition', refit_cooks: bool=True) -> pd.DataFrame:
	counts, metadata = metacounts
	print(f"\t\tRunning DDS on {strain}...")
	dds = _dds_strain(counts, metadata, design_factors, refit_cooks)
	print(f"\t\tRunning DS on {strain}...")
	degs = _ds_strain(dds, strain)
	print(f"\t\tFinished processing {strain}.")

	return degs


def _mprundeseq(strains: list, metacounts: list[tuple[pd.DataFrame, pd.DataFrame]], nworkers: int) -> list[pd.DataFrame]: # idc about dds params rn
	exec = ProcessPoolExecutor(max_workers=nworkers)
	degs = list(exec.map(_rundeseq, strains, metacounts))
	exec.shutdown(wait=False, cancel_futures=True)
	# degs = []
	# for strain, metacount in zip(strains, metacounts):
	# 	degs.append(_rundeseq(strain, metacount))

	return degs

def _evalsignificance_strain(deg: pd.DataFrame, l2fc: float, padj: float) -> dict[str: int]:
	sigmap = {}

	for i, row in deg.iterrows():
		transcript = i
		if row['padj'] < padj and abs(row['log2FoldChange']) > l2fc: # significant
			if row['log2FoldChange'] < 0:
				sigmap[transcript] = -1
			elif row['log2FoldChange'] > 0:
				sigmap[transcript] = 1
		else:
			sigmap[transcript] = 0

	return sigmap


def _deseq_pipeline(tpms: pd.DataFrame, strains: list, l2fc: float, padj: float, nworkers: int, export: Path) -> pd.DataFrame:
	sigmatrix = pd.DataFrame(columns=['transcript', 'N2'] + strains)
	sigmatrix['transcript'] = tpms.index
	sigmatrix['N2'] = sigmatrix['N2'].fillna(0)
	
	print("\tGetting strain-wise metadata and counts...")
	metacounts = [_metacount(tpms, strain) for strain in strains]

	print("\tRunning DeSeq2 analysis...")
	start = time.perf_counter()
	degs = _mprundeseq(strains, metacounts, nworkers=nworkers)

	end = time.perf_counter()
	print(f"\t{(end - start):.6f} seconds elapsed.")

	if export:
		for strain, deg in zip(strains, degs):
			name = f"{strain}_N2_DEGs"
			deg.to_csv(export/f'{name}.tsv', sep='\t')

	print("\tEvaluating significance...")
	for strain, deg in zip(strains, degs):
		sigmap = _evalsignificance_strain(deg, l2fc, padj)
		sigmatrix[strain] = sigmatrix['transcript'].map(sigmap)

	return sigmatrix


def sigmatrix(tpmpath: Path, loss_strains: list=[], keep_strains: list=[], l2fc: float=1.5, padj: float=0.05, nworkers: int=4, export: Path=None) -> pd.DataFrame:
	print("Reading TPMs...")
	tpms = pd.read_csv(tpmpath, sep='\t')

	print("Filtering for strains...")
	tpms = _targetstrains(tpms, loss_strains, keep_strains)

	print("Running DeSeq2 pipeline...")
	strains = loss_strains + keep_strains
	sigmatrix = _deseq_pipeline(tpms, strains, l2fc, padj, nworkers, export=export)
	print("Done!")
	return sigmatrix


def _setop(row: pd.Series, loss_strains: list, keep_strains: list) -> bool:
	loss_set = set()
	keep_set = set()
	for strain in loss_strains:
		loss_set.add(row[strain])

	for strain in keep_strains:
		keep_set.add(row[strain])

	if 0 not in loss_set and not loss_set.intersection(keep_set):
		return True
	else:
		return False


def findsigdegs(degdir: Path, sigmatrixpath: Path, loss_strains: list, keep_strains: list, export: Path=None) -> list:
	# tpms = pd.read_csv(tpmpath, sep='\t').reset_index(drop=True)
	degs = {}
	for degfile in degdir.iterdir():
		if degfile.is_file():
			strain = degfile.stem.split('_')[0]
			degs[strain] = pd.read_csv(degfile, sep='\t')
	
	sigmatrix = pd.read_csv(sigmatrixpath, sep='\t')

	sigtranscripts = []

	for i, row in sigmatrix.iterrows():
		op = _setop(row, loss_strains, keep_strains)
		if op:
			sigtranscripts.append(row['transcript'])
	
	if export:
		for strain, deg in degs.items():
			name = f"{strain}_candidates"
			straindf = deg.loc[deg['transcript'].isin(sigtranscripts)]
			straindf.to_csv(export/f'{name}.tsv', sep='\t', index=False)
	
	return sigtranscripts


if __name__ == "__main__":
	DIR = Path("/Users/kvying/Files/ucsb-local/rothman/deg-pipeline")
	OUTDIR = DIR / "out"
	DEGDIR = OUTDIR / 'degs'
	CANDIDATESDIR = OUTDIR / 'candidates'
	os.makedirs(OUTDIR, exist_ok=True)
	os.makedirs(DEGDIR, exist_ok=True)
	os.makedirs(CANDIDATESDIR, exist_ok=True)

	loss_strains = ['CX11314']
	keep_strains = ['LKC34', 'ED3017', 'JU775', 'MY16', 'MY23', 'JT11398', 'CB4856']
	sm = sigmatrix(OUTDIR/'resolved_tpms.tsv', loss_strains, keep_strains, nworkers=4, export=DEGDIR)
	sm.to_csv(OUTDIR/'sigmatrix.tsv', sep='\t', index=False)

	sigdegs = findsigdegs(DEGDIR, OUTDIR/'sigmatrix.tsv', loss_strains, keep_strains, CANDIDATESDIR)
	
