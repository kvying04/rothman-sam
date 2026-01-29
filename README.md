# Sam Stuff

## `deg-pipeline`
**`isoforms.py`**: Resolves discrepancies in naming of transcripts in RNAseq data by cross-referencing each isoform with Wormbase through its REST API. 

**`intersections.py`**: Finds significant DEGs that serve as candidate genes for an observed phenotype. First finds differentially expressed genes (DEGs) and thresholds significance using DESeq2. From this, generates a significance matrix encoding directionality of log-fold change from N2 (if $p_{adj} < \alpha$). Finds genes in significance matrix consistently down/upregulated in the `loss` group but variant in the `keep` group. 

## `snps-pipeline`
**`getsnps.py`**:

**`snpstogene.py`**:

**`geneoverlaps.py`**:

## `genepaircorr`
**`spearmancorr.py`**: Stitches TPMs generated from `aligngtf.sh` into single Pandas dataframe. Calculates a Spearman rank correlation coefficient for each gene pair. P-values are adjusted for false discovery rate (FDR) using Benjamini-Hochberg and Benjamini-Yekutieli procedures. 
### `0-srx-to-srr`
**`srxclean.py`**: Parses gene expression accession ID table from Wormbase into a `.txt` file. 

**`getsrrlists.sh`**: Queries using the SRR lists of all provided SRX's using Entrez Direct. 

**`stitchsrrlists.sh`**: Concatenates all SRR lists produced from `getsrrlists.sh` into one `.txt` file. 
### `1-srr-to-fastq`
**`prefetch.sh`**:

**`getfastqs.sh`**:
### `2-fastq-to-bam`
**`alignfasta.sh`**:
### `3-bam-to-tpm`
**`aligngtf.sh`**: Aligns BAMs produced from `alignfasta.sh` against a GTF annotation file (in this case we used WBcel235/ce11 reference genome from NCBI). Produces `.tpm` files.