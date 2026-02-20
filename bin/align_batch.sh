#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = fasta_files
# $4 = file_type

touch alignment.cfa
for FILE in ${3}; do
    # do alignment only if number of sequences is >1 
    if [[ $( grep -c "^>" $FILE ) > 1 ]]; then
        if [[ $4 == "small" ]]; then
            # align with highly accurate, slower mode for small files
            mafft \
                --nuc \
                --thread ${2} \
                --linelength -1 \
                --globalpair \
                --maxiterate 1000 \
                $FILE | \
            sed -e ':a;N;$!ba;s/\n/>>>/g' \
            >> alignment.cfa
        else 
            # align with faster, slightly less accurate method for larger files
            mafft \
                --nuc \
                --thread ${2} \
                --linelength -1 \
                --maxiterate 10 \
                $FILE | \
            sed -e ':a;N;$!ba;s/\n/>>>/g' \
            >> alignment.cfa
        fi
    elif [[ $( grep -c "^>" $FILE ) == 1 ]]; then
        # else rename input file as output file (no alignment needed)
        echo "No alignment required for $FILE"
        cat $FILE |
        sed -e ':a;N;$!ba;s/\n/>>>/g' \
        >> alignment.cfa 
    else 
        # else don't output anything (empty file)
        echo "Empty file input -- not outputting a file"
    fi
done
