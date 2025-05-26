#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = translations
# $4 = hmm

# replace each space in sequence headers with the string "!?!?"
sed '/^>/ s/ /!?!?/g' $3 > renamed.fasta

# find hits
hmmsearch \
    --domtblout hmmer_domtblout.txt \
    --notextw \
    --cpu $2 \
    $4 \
    renamed.fasta \
    > hmmer.out

### TODO: add parameters for search