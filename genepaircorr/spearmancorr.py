import pandas as pd
from scipy.stats import spearmanr
from pathlib import Path
from statsmodels.stats.multitest import multipletests

ROOT = Path("/home/hanwenying/rothman/genepaircorr")
TPMDIR = ROOT/"3-bam-to-tpm/tpms"

ALPHA = 0.05

# stitch tpms
for i, tpmpath in enumerate(TPMDIR.glob('*.tpm')):
	srrcol = pd.read_csv(tpmpath, sep='\t') # DOUBLE CHECK DELMITER
	if i == 0:
		tpms = srrcol
	else:
		srrid = tpmpath.stem
		tpms[srrid] = srrcol

# spearman correlation: genes x srr shape
rhos, pvals = spearmanr(tpms.T) # 2 dataframes of g x g shape

pmatrix = pvals.to_numpy()
pshape, pvec = pmatrix.shape, pmatrix.flatten()

# FDR correction: use Benjamini-Hochberg AND Benjamini-Yekutieli
padj_bh = multipletests(pvec, alpha=ALPHA, method='fdr_bh', is_sorted=False, returnsorted=False)
padj_by = multipletests(pvec, alpha=ALPHA, method='fdr_by', is_sorted=False, returnsorted=False)

pmatrix_bh = padj_bh[1].reshape(pshape)
pmatrix_by = padj_by[1].reshape(pshape)

