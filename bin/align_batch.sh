#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = fasta_files
# $4 = file_type


for FILE in ${3}; do
    # get basename of file
    BASE=$( basename "$FILE" .fasta )
    # do alignment only if number of sequences is >1 
    if [[ $( grep -c "^>" $FILE ) > 1 ]]; then
        if [[ $4 == "small" ]]; then
            # align with highly accurate, slower mode for small files
            mafft \
                --nuc \
                --thread ${2} \
                --linelength -1 \
                --globalpair \
                --maxiterate 10 \
                $FILE \
                > ${BASE}.aligned.fasta
        else 
            # align with faster, slightly less accurate method for larger files
            mafft \
                --nuc \
                --thread ${2} \
                --linelength -1 \
                --maxiterate 10 \
                $FILE \
                > ${BASE}.aligned.fasta
        fi
    elif [[ $( grep -c "^>" $FILE ) == 1 ]]; then
        # else rename input file as output file (no alignment needed)
        echo "No alignment required for $FILE"
        cp $FILE ${BASE}.aligned.fasta
    else 
        # else don't output anything (empty file)
        echo "Empty file input -- not outputting a file"
    fi
done
