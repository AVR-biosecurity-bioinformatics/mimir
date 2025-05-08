#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = fasta_file

# remove exact duplicates based on name
seqkit rmdup \
    --by-name \
    --line-width 0 \
    --threads $2 \
    --dup-seqs-file removed.fasta \
    < $3 \
    > seqs_deduplicated.fasta

# touch removed.fasta if doesn't exist
if [ -f removed.fasta ]; then
    echo "Finished deduplicating sequences"
else 
    touch removed.fasta
    echo "Finished deduplicating sequences -- 'removed.fasta' created directly"
fi 