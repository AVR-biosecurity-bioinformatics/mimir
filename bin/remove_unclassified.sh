#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = fasta_file
# $4 = '--remove_unclassified' parameter
# $5 = task.index

OUTPUT_FILE="remove_unclassified.${5}.fasta"

case ${4} in 

	"all_ranks")
		seqkit grep -n -r -p "(;Unclassified){7}$" -v -w 0 \
	    $3 \
	    -o $OUTPUT_FILE
		;;
	
	"any_ranks")
		seqkit grep -n -r -p ";(.*)?Unclassified(.*)?$" -v -w 0 \
	    $3 \
	    -o $OUTPUT_FILE
		;;

    "terminal")
        seqkit grep -n -r -p "Unclassified$" -v -w 0 \
	    $3 \
	    -o $OUTPUT_FILE
        ;;
	
	"none")
	    cp $3 $OUTPUT_FILE
		;;
	
	*)
		echo "'--remove_unclassified' value $REMOVE_UNCLASSIFIED not valid "
		exit 1
		;;
esac
