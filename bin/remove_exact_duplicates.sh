#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = fasta_file

# remove exact duplocates
seqkit rmdup \
    --by-name \
    --line-width 0 \
    < $3 \
    > seqs_deduplicated.fasta

