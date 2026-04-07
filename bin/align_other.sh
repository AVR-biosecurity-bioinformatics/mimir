#!/bin/bash
set -e
set -u
set -o pipefail
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = core.fasta
# $4 = other.fasta

# if 'other' sequence file is empty, skip alignment
if [ -s $4 ]; then 
    # align
    mafft \
        --nuc \
        --thread ${2} \
        --linelength -1 \
        --add ${4} \
        ${3} | \
    sed -e ':a;N;$!ba;s/\n/>>>/g' \
    > alignment.cfa 
else 
    # process input as .cfa format
    echo "All sequences are already aligned"
    cat $3 |
    sed -e ':a;N;$!ba;s/\n/>>>/g' \
    > alignment.cfa 
fi 

# throw error if output file is empty
if [ -s alignment.cfa  ]; then
    echo "Finished aligning all sequences"        
else 
    echo "alignment output file is empty"
    exit 1
fi
