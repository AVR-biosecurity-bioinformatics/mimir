#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = translations
# $4 = hmm


# find hits
hmmsearch \
    --domtblout hmmer_domtblout.txt \
    --notextw \
    --cpu $2 \
    $4 \
    $3 \
    > hmmer.out

### TODO: add parameters for search