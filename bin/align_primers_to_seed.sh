#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = fasta_file
# $4 = primers


# align primers to database alignment
mafft \
    --amino \
    --thread $2 \
    --linelength -1 \
    --keeplength \
    --addfragments $4 \
    "${3}" \
    > aligned_primers.fasta

# throw error if output file is empty
if [ -s aligned_primers.fasta ]; then
    echo "Finished aligning primer sequences"        
else 
    echo "alignment output file is empty"
    exit 1
fi