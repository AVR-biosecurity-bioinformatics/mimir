#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = fasta_file

# concat chunked .fasta files together
# NOTE: Needs all files to be staged with "path()" process input
cat *.fasta > chunks_combined.fasta

# dealign sequences
seqkit seq \
    --gap-letters "-" \
    --remove-gaps \
    --line-width 0 \
    < chunks_combined.fasta \
    > chunks_combined_dealigned.fasta

