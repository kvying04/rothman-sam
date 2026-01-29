ROOT="/home/hanwenying/rothman/genepaircorr"
SRXDIR="$ROOT/0-srx-to-srr"
SRRDIR="$ROOT/1-srr-to-fastq"
PREFETCHDIR="$SRRDIR/prefetch"
srrlist="$SRXDIR/totalsrrs.txt"

while IFS= read -r line; do
	echo $line
    prefetch $line -O $PREFETCHDIR
done < "$srrlist"