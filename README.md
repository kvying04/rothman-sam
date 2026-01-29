# Sam Stuff

## `deg-pipeline`
**`isoforms.py`**: Resolves discrepancies in naming of transcripts in RNAseq data by cross-referencing each isoform with Wormbase through its REST API. 

**`intersections.py`**: Finds significant DEGs that serve as candidate genes for an observed phenotype. First finds differentially expressed genes (DEGs) and thresholds significance using PyDESeq2.[^1] From this, generates a significance matrix encoding directionality of log-fold change from N2 (if $p_{adj} < \alpha$). Finds genes in significance matrix consistently down/upregulated in the loss group but variant in the keep group. 

## `snps-pipeline`
**`getsnps.py`**:

**`snpstogene.py`**:

**`geneoverlaps.py`**:

## `genepaircorr`
**`spearmancorr.py`**: Stitches TPMs generated from `aligngtf.sh` into single Pandas dataframe. Calculates a Spearman rank correlation coefficient for each gene pair.[^2] P-values are adjusted for false discovery rate (FDR) using Benjamini-Hochberg and Benjamini-Yekutieli procedures.[^3]
### `0-srx-to-srr`
**`srxclean.py`**: Parses gene expression accession ID table[^4] from Wormbase into a text file. 

**`getsrrlists.sh`**: Queries using the SRR lists of all provided SRX's using Entrez Direct.[^5]

**`stitchsrrlists.sh`**: Concatenates all SRR lists produced from `getsrrlists.sh` into one text file. 
### `1-srr-to-fastq`
**`prefetch.sh`**:

**`getfastqs.sh`**:
### `2-fastq-to-bam`
**`alignfasta.sh`**:
### `3-bam-to-tpm`
**`aligngtf.sh`**: Aligns BAMs produced from `alignfasta.sh` against a GTF annotation file (in this case we used the WBcel235/ce11 reference genome[^6] from NCBI). Produces TPM files.


[^1]: [PyDESeq2: a python package for bulk RNA-seq differential expression analysis.](https://pydeseq2.readthedocs.io/en/stable/)
[^2]: [spearmanr — SciPy v1.17.0 Manual.](https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.stats.spearmanr.html)
[^3]: [statsmodels.stats.multitest.multipletests](https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html)
[^4]: No idea how to find a proper accession to it...
[^5]: Kans J. Entrez® Direct: E-utilities on the Unix Command Line. 2013 Apr 23 [Updated 2025 Mar 25]. In: Entrez® Programming Utilities Help [Internet]. Bethesda (MD): National Center for Biotechnology Information (US); 2010-. Available from: https://www.ncbi.nlm.nih.gov/books/NBK179288/
[^6]: [Genome assembly WBcel235.](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002985.6/)