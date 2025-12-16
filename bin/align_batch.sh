#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = fasta_files

# $3 is a list of .fasta files

# # convert groovy list to bash list
# FILE_LIST=$(echo $3 | tr -d '[],')

for FILE in ${3}; do
    # get basename of file
    BASE=$( basename "$FILE" .fasta )
    # do alignment only if number of sequences is >1 
    if [[ $( grep -c "^>" $FILE ) > 1 ]]; then
        # align
        mafft \
            --nuc \
            --thread ${2} \
            --linelength -1 \
            --globalpair \
            --maxiterate 10 \
            $FILE \
            > ${BASE}.aligned.fasta
    else 
        # else rename input file as output file (no alignment needed)
        echo "No alignment required for $FILE"
        cp $FILE ${BASE}.aligned.fasta
    fi
done
