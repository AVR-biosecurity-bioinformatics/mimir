#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = fasta_file

# $3 is a list of .fasta files

for FILE in $3
do
    # get basename of file
    BASE=$( basename "$FILE" .fasta )
    # do alignment only if number of sequences is >1 AND sequence headers don't contain ";Unclassified"
    if [[ $( grep -c "^>" $FILE ) > 1 ]] && [[ $( head $FILE -n1 | grep ";Unclassified" | wc -l ) = 0 ]]; then
        clustalo \
            -i $FILE \
            -o ${BASE}.aligned.fasta \
            --outfmt=fasta \
            --wrap=999999 \
            --seqtype=DNA \
            --dealign \
            --threads=${2}
    else 
        # else rename input file as output file (no alignment needed)
        cp $FILE ${BASE}.aligned.fasta
    fi
done


