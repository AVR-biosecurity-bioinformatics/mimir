#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = fasta_file
# $4 = subMSAtable

### for family-level alignment
export MAFFT_TMPDIR="/tmp"

# if only one sequence in file, skip alignment
if [[ $( grep -c "^>" $3 ) > 1 ]]; then 
    # align
    mafft \
        --nuc \
        --thread $2 \
        --linelength -1 \
        --merge $4 \
        "${3}" \
        > aligned.fasta
else 
    # rename as aligned
    cp $3 aligned.fasta
fi 
