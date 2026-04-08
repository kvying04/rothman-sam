ROOT="/home/hanwenying/rothman-sam/w76"
OUTDIR="$ROOT/out"
VCFDIR="$OUTDIR/vcfs"
ZIPDIR="$OUTDIR/vcfsgz"

for vcfpath in $VCFDIR/*.vcf; do
	sample=$(basename "$vcfpath" .vcf)
	bgzip $vcfpath -o $ZIPDIR/$sample.vcf.gz 
	bcftools index $ZIPDIR/$sample.vcf.gz
done

ls $ZIPDIR/*.vcf.gz > $VCFDIR/vcflist.txt

bcftools merge -l $VCFDIR/vcflist.txt -Ov -o $OUTDIR/merged.vcf