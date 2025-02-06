#!/bin/bash
set -e
set -u
set -o pipefail
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

# throw error if output file is empty
if [ -s aligned.fasta ]; then
    echo "Finished aligning core sequences"        
else 
    echo "alignment output file is empty"
    exit 1
fi