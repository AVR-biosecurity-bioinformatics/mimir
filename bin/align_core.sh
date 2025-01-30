#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = core.fasta

# if only one sequence in file, skip alignment
if [[ $( grep -c "^>" $3 ) > 1 ]]; then 
    # align
    mafft \
        --nuc \
        --thread ${2} \
        --linelength -1 \
        --globalpair \
        --maxiterate 10 \
        ${3} \
        > aligned.fasta
else 
    # rename as aligned
    cp $3 aligned.fasta
fi 
