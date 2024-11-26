#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = fasta_file
# $4 = split_rank

## parse args

# check input fasta is unwrapped, and if wrapped, unwrap
if head -n3 $3 | grep "^>" | wc -l | grep -q 2; then
	# if first three lines contain two headers, input is unwrapped
	INPUT_FASTA=${3}
else 
	# if first three lines don't contain two, input is likely wrapped
	# unwrap fasta
	awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $3 > unwrapped.fasta
	INPUT_FASTA=unwrapped.fasta
fi

# get unique lineage strings for the split rank
SPLIT_RANK="${4}"

case $SPLIT_RANK in 
	"species")
		REGEX=0
		;;
	"genus")
		REGEX=1
		;;
	"family")
		REGEX=2
		;;
	"order")
		REGEX=3
		;;
	"class")
		REGEX=4
		;;
	"phylum")
		REGEX=5
		;;
	"kingdom")
		REGEX=6
		;;
	*)
		echo "SPLIT_RANK '$SPLIT_RANK' not valid"
		exit 1
		;;
esac

grep -Po "(?<=;).+?(?=(;[^;]+){$REGEX}$)" $INPUT_FASTA | sort -u > lineage_strings.txt

# loop through lineage strings, saving all matching sequences to a file
while IFS="" read -r p || [ -n "$p" ]
do
	grep \
	--no-group-separator \
	-A1 \
	-E ";${p}" \
	$INPUT_FASTA \
	> ${p//;/_}.lineage.fasta
done < lineage_strings.txt

