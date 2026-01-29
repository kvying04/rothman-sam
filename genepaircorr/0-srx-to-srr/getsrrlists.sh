ROOT="/home/hanwenying/rothman/genepaircorr"
SRXDIR="$ROOT/0-srx-to-srr"
SRRDIR="$SRXDIR/srrlists"

mkdir -p "$SRRDIR"

srxs="$SRXDIR/srxlist.txt"
TOTAL=$(wc -l < "$srxs")

while read -r line; do
	[[ -z "$line" ]] && continue
    esearch -db sra -query "$line" < /dev/null | efetch -format runinfo | cut -d ',' -f 1 | grep -E '^[S|E|D]RR' > "${SRRDIR}/${line}_SRRlist.txt"
	echo "${line}"
done < "$srxs" | pv -l -s $TOTAL > /dev/null