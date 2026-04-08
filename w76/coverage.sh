DATA="/home/hanwenying/rothman-data/w76/bams"
OUT="/home/hanwenying/rothman-sam/w76/out"

for bampath in $DATA/*.bam; do
    sample=$(basename "$bampath" .bam)
    samtools coverage "$bampath" > "${OUT}/coverages/${sample}-coverage.tsv"
done