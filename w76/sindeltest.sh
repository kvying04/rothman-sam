DATA="/home/hanwenying/rothman-data/w76/bams"
OUT="/home/hanwenying/rothman-sam/w76/imsindel-out"
FASTAPATH="/home/hanwenying/rothman-sam/ref/ce11_genome.fna"

for FILE in $DATA/*.bam; do
	for chr in {"NC_003279.8","NC_003280.10","NC_003281.10","NC_003282.8","NC_003283.11","NC_003284.9","NC_001328.1"}; do
		# docker run --rm -v $DATA:/data 
		imsindel --bam $FILE --chr $chr --outd "$OUT/$(basename $FILE)" --indelsize 10000 --reffa $FASTAPATH --thread 8
	done
done

# FILE=$DATA/I9.bam
# for chr in {"NC_003279.8","NC_003280.10","NC_003281.10","NC_003282.8","NC_003283.11","NC_003284.9","NC_001328.1"}; do
# 	imsindel --bam $FILE --chr $chr --outd "$OUT/$(basename $FILE)" --indelsize 10000 --reffa $FASTAPATH --thread 8
# done

# FILE=$DATA/I9.bam
# chr="NC_001328.1"
# imsindel --bam $FILE --chr $chr --outd "$OUT/$(basename $FILE)" --indelsize 10000 --reffa $FASTAPATH --thread 8