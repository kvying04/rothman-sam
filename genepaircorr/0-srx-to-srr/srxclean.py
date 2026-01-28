import pandas as pd
from pathlib import Path

DIR = Path("/Users/kvying/Files/ucsb-local/rothman/genepaircorr")
SRXDIR = DIR/"0-srx-to-srr"
OUTPATH = SRXDIR/"srxlist.txt"

rawsrx = pd.read_csv(SRXDIR/"wbsrx.tsv", sep='\t')
srxs = list(rawsrx['Name'].apply(lambda name: name.split('.')[-1]))

with open(OUTPATH, 'w') as f:
	for srx in srxs:
		f.write(f"{srx}\n")

print(f"SRX list saved to {OUTPATH}")