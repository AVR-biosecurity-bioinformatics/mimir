#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = fasta_file

clustalo \
    -i $3 \
    -o aligned.fasta \
    --outfmt=fasta \
    --wrap=0 \
    --seqtype=DNA \
    --dealign \
    --threads=${2}
