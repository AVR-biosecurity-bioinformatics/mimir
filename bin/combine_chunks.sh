#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = fasta_file

# concat chunked .fasta files together
# NOTE: Needs all files to be staged with "path()" process input
cat *.fasta > chunks_combined_pre.fasta

if [ "$4" == "true" ]; then
    # dealign sequences
    seqkit seq \
        --gap-letters "-" \
        --remove-gaps \
        --line-width 0 \
        --threads $2 \
        < chunks_combined_pre.fasta \
        > chunks_combined.fasta
else 
    mv chunks_combined_pre.fasta chunks_combined.fasta
fi

