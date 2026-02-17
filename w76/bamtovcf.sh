DATA="/home/hanwenying/rothman-data/w76/bams"
OUT="/home/hanwenying/rothman-sam/w76/imsindel-out"
FASTAPATH="/home/hanwenying/rothman-sam/ref/ce11_genome.fna"

for bampath in $DATA/*.bam; do
	sample=$(basename "$bampath" .bam)
	echo $sample
	bcftools mpileup -Ou -f $FASTAPATH "$bampath" | bcftools call -mv -Ob -o "${sample}".bcf
	bcftools view "${sample}".bcf -o "${sample}".vcf
	bcftools index "${sample}".vcf
done