#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = fasta_file
# $4 = split_rank
# $5 = task.time (dynamic, in seconds)

## parse args
N_SEQS=$( grep -c "^>" $3 )

echo "Dynamic task time for $N_SEQS sequences is ${5} seconds"
# set seconds to 0 
SECONDS=0

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

# get unique lineage strings for the split rank
SPLIT_RANK="${4}"

### new implementation with hashed filenames

# set the number of ; separators to join together when defining the lineage string
case $SPLIT_RANK in 
	"kingdom")
		N_SEPS="2"
		;;
	"phylum")
		N_SEPS="3"
		;;
	"class")
		N_SEPS="4"
		;;
	"order")
		N_SEPS="5"
		;;
	"family")
		N_SEPS="6"
		;;
	"genus")
		N_SEPS="7"
		;;
	"species")
		N_SEPS="8"
		;;	
	*)
		echo "Invalid option for '$SPLIT_RANK': $SPLIT_RANK"
		exit 1
		;;
esac

# save sequences to new files based on their lineage string hash
awk -v n_seps="$N_SEPS" -v RS=">" -v ORS="" -v FS="[\r\n]+" -v OFS="\n" ' NR>1 {
    split($1, he, ";", seps )
    lineage_string=""
    for (i = 2; i <= n_seps; i++ )
        lineage_string = lineage_string ";" he[i]
    tmp="echo \""lineage_string"\" | md5sum | cut -f1 -d\" \""
    tmp | getline linhash
    gsub(/[^[:alnum:]]/, "", linhash)
    print ">" $0 > ( linhash ".lineage.fasta")
    }' \
    < $FASTA

# compare elapsed and requested time 
echo "${SECONDS} seconds elapsed; ${5} requested."

### old implementation
# case $SPLIT_RANK in 

# 	"kingdom")
# 		awk -v RS=">" -v ORS="" -v FS="[\r\n]+" -v OFS="\n" 'NR>1 { split($1, he, ";", seps ); print ">" $0 > (he[2] ".lineage.fasta") }' $FASTA
# 		;;
# 	"phylum")
# 		awk -v RS=">" -v ORS="" -v FS="[\r\n]+" -v OFS="\n" 'NR>1 { split($1, he, ";", seps ); print ">" $0 > (he[2] "_" he[3] ".lineage.fasta") }' $FASTA
# 		;;
# 	"class")
# 		awk -v RS=">" -v ORS="" -v FS="[\r\n]+" -v OFS="\n" 'NR>1 { split($1, he, ";", seps ); print ">" $0 > (he[2] "_" he[3] "_" he[4] ".lineage.fasta") }' $FASTA
# 		;;
# 	"order")
# 		awk -v RS=">" -v ORS="" -v FS="[\r\n]+" -v OFS="\n" 'NR>1 { split($1, he, ";", seps ); print ">" $0 > (he[2] "_" he[3] "_" he[4] "_" he[5] ".lineage.fasta") }' $FASTA
# 		;;
# 	"family")
# 		awk -v RS=">" -v ORS="" -v FS="[\r\n]+" -v OFS="\n" 'NR>1 { split($1, he, ";", seps ); print ">" $0 > (he[2] "_" he[3] "_" he[4] "_" he[5] "_" he[6] ".lineage.fasta") }' $FASTA
# 		;;
# 	"genus")
# 		awk -v RS=">" -v ORS="" -v FS="[\r\n]+" -v OFS="\n" 'NR>1 { split($1, he, ";", seps ); print ">" $0 > (he[2] "_" he[3] "_" he[4] "_" he[5] "_" he[6] "_" he[7] ".lineage.fasta") }' $FASTA
# 		;;
# 	"species")
# 		awk -v RS=">" -v ORS="" -v FS="[\r\n]+" -v OFS="\n" 'NR>1 { split($1, he, ";", seps ); print ">" $0 > (he[2] "_" he[3] "_" he[4] "_" he[5] "_" he[6] "_" he[7] "_" he[8] ".lineage.fasta") }' $FASTA
# 		;;	
# 	*)
# 		echo "Invalid option for '$SPLIT_RANK': $SPLIT_RANK"
# 		exit 1
# 		;;
# esac

# # compared elapsed and requested time 
# echo "${SECONDS} seconds elapsed; ${5} requested."