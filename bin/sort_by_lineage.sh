#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = fasta_file
# $4 = task.time (dynamic, in seconds)

## parse args
N_SEQS=$( grep -c "^>" $3 )

echo "Dynamic task time for $N_SEQS sequences is ${4} seconds"
# set seconds to 0 
SECONDS=0

# check input fasta is unwrapped, and if wrapped, unwrap
if head -n3 $3 | grep "^>" | wc -l | grep -q 2; then
	# if first three lines contain two headers, input is unwrapped
	FASTA=${3}
else 
	# if first three lines don't contain two, input is likely wrapped
	# unwrap fasta
	awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $3 > unwrapped.fasta
	FASTA=unwrapped.fasta
fi

### sort file by lineage string
## inspired by https://unix.stackexchange.com/questions/602420/sort-fasta-file-based-on-its-alphanumeric-ids
# remove null characters if present
# then add lineage string to start of header with seven ; as a spacer
# then add null character between records
# then sort records 
# then remove null characters
# then remove lineage string from the start of header

tr < $FASTA -d '\000' \
    | awk -v RS=">" -v ORS="\x00" -v FS="[\r\n]+" -v OFS="\n" 'NR>1 { split($1, he, ";", seps ); print ">" he[2] "_" he[3] "_" he[4] "_" he[5] "_" he[6] "_" he[7] "_" he[8] "|||||||" $0 }' - \
    | sort -z \
    | tr -d '\000' \
    | sed --expression=$'s/^>.*|||||||/>/'\
    > sorted.fasta

# compared elapsed and requested time 
echo "${SECONDS} seconds elapsed; ${4} requested."