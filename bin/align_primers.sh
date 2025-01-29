#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = fasta_file
# $4 = primers


##


# align primers to database alignment
mafft \
    --nuc \
    --thread $2 \
    --linelength -1 \
    --6merpair \
    --keeplength \
    --addfragments $4 \
    "${3}" \
    > aligned_primers.fasta
