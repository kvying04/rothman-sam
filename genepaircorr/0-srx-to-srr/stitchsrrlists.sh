ROOT="/home/hanwenying/rothman/genepaircorr"
SRXDIR="$ROOT/0-srx-to-srr"
SRRDIR="$SRXDIR/srrlists"

# # while read -r line; do
# # 	[[ -z "$line" ]] && continue
# #     esearch -db sra -query "$line" < /dev/null | efetch -format runinfo | cut -d ',' -f 1 | grep SRR > "${SRRDIR}/${line}_SRRlist.txt"
# # 	echo "${line}"
# # done < "$srxs" | pv -l -s $TOTAL > /dev/null

# for file in $SRRDIR; do

cat "$SRRDIR"/[SED]*.txt >> totalsrrs.txt

# cat "$SRR"/[SED]*.txt >> totalsrrs.txt