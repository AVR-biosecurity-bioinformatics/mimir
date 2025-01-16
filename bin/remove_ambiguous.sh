#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = fasta_file

# check input fasta is unwrapped, and if wrapped, unwrap
if head -n3 $3 | grep "^>" | wc -l | grep -q 2; then
	# if first three lines contain two headers, input is unwrapped
	FASTA=${3}
else 
	# if first three lines don't contain two headers, input is likely wrapped
	# unwrap fasta
	awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $3 > unwrapped.fasta
	FASTA=unwrapped.fasta
fi

# remove all sequences that contain any bases other than ACGT
cat $FASTA | \
	awk -vRS=">" -vORS="" -vFS="\n" -vOFS="\n" '
    	NR>1 && $2!~/[RYSWKMBDHVNIUryswkmbdhvniu]/ {print ">"$0}
  		' - \
	> remove_ambiguous.fasta

# capture removed sequences
cat $FASTA | \
	awk -vRS=">" -vORS="" -vFS="\n" -vOFS="\n" '
    	NR>1 && $2!~/[RYSWKMBDHVNIUryswkmbdhvniu]/ {print ">"$0}
  		' - \
	> removed.fasta
