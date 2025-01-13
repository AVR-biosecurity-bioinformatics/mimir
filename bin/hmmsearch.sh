#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads


# find hits
hmmsearch \
    --domtblout hmmer_domtblout.txt \
    --notextw \
    hmm.gz \
    translations.fasta \
    > hmmer.out

### TODO: add parameters for search