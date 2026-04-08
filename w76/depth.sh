DATA="/home/hanwenying/rothman-data/w76/bams"
OUT="/home/hanwenying/rothman-sam/w76/out"


for bampath in $DATA/*.bam; do
	sample=$(basename "$bampath" .bam)
	samtools depth -a "$bampath" > "${OUT}/depths/${sample}-all.tsv"
	# samtools depth -a -r NC_003280.10:5833328-5836203 "$bampath" > "${OUT}/depths/${sample}-fzo1.tsv"
done