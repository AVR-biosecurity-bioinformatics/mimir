#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = core.aligned.fasta
# $4 = other.fasta


# if 'other' sequence file does not contain sequences, skip alignment
if [[ $( grep -c "^>" $4 ) > 0 ]]; then 
    # align
    mafft \
        --nuc \
        --thread ${2} \
        --linelength -1 \
        --memsavetree \
        --add ${4} \
        ${3} \
        > all.aligned.fasta
else 
    # rename as aligned
    cp $3 all.aligned.fasta
fi 

# throw error if output file is empty
if [ -s all.aligned.fasta ]; then
    echo "Finished aligning all sequences"        
else 
    echo "alignment output file is empty"
    exit 1
fi