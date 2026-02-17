DATA="/home/hanwenying/rothman-data/w76/bams"
OUT="/home/hanwenying/rothman-sam/w76/imsindel-out"
FASTAPATH="/home/hanwenying/rothman-sam/ref/ce11_genome.fna"

for FILE in $DATA/*.bam; do
	for chr in {'I', 'II', 'III', 'IV', 'V', 'X', 'MtDNA'}; do
		# docker run --rm -v $DATA:/data 
		imsindel --bam $FILE --chr $chr --outd $OUT --indelsize 10000 --reffa $FASTAPATH
	done
done