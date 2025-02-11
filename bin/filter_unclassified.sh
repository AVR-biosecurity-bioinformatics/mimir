#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = fasta_file
# $4 = '--remove_unclassified' parameter

OUTPUT_FILE="filter_unclassified.fasta"
REMOVED_FILE="removed.fasta"

case ${4} in 

	"all_ranks")
		# retained sequences
		seqkit grep -n -r -p "(;Unclassified){7}$" -v -w 0 \
	    	$3 \
	    	-o $OUTPUT_FILE
		
		# removed sequences
		seqkit grep -n -r -p "(;Unclassified){7}$" -w 0 \
	    	$3 \
	    	-o $REMOVED_FILE
		;;
	
	"any_ranks")
		# retained sequences
		seqkit grep -n -r -p ";(.*)?Unclassified(.*)?$" -v -w 0 \
	    	$3 \
			-o $OUTPUT_FILE
		
		# removed sequences
		seqkit grep -n -r -p ";(.*)?Unclassified(.*)?$" -w 0 \
	    	$3 \
	    	-o $REMOVED_FILE
		;;

    "terminal")
		# retained sequences
        seqkit grep -n -r -p "Unclassified$" -v -w 0 \
	    	$3 \
	    	-o $OUTPUT_FILE
        
		# removed sequences
		seqkit grep -n -r -p "Unclassified$" -w 0 \
	    	$3 \
	    	-o $REMOVED_FILE
		;;
	
	"none")
		# retained sequences
	    cp $3 $OUTPUT_FILE

		# removed sequences
		touch $REMOVED_FILE
		;;
	
	*)
		echo "'--remove_unclassified' value $REMOVE_UNCLASSIFIED not valid "
		exit 1
		;;
esac
