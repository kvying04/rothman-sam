import pandas as pd

loss_strains = ['CX11314']
keep_strains = ['LKC34', 'ED3017', 'JU775', 'MY16', 'MY23', 'JT11398', 'CB4856']

isolates = loss_strains + keep_strains

for isolate in isolates:
	inpath = f"/Users/kvying/Files/ucsb-local/rothman-sam/deg-pipeline/out/degs/{isolate}_N2_DEGs.tsv"

	degs = pd.read_csv(inpath, sep='\t')
	degs = degs[degs['baseMean'] != 0]
	degs = degs[abs(degs['log2FoldChange']) > 1.5]
	degs = degs[degs['padj'] < 0.05]

	degs.to_csv(f"/Users/kvying/Files/ucsb-local/rothman-sam/deg-pipeline/out/shinygo/forkegg_{isolate}.tsv", index=False, sep='\t')