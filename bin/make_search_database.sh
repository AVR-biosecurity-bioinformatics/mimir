#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = target_fasta

# replace each space in sequence headers of fasta file with the string "!?!?"
sed '/^>/ s/ /!?!?/g' $3 > target.fasta

# make database
lastdb \
    -v \
    -P $2 \
    -S 1 \
    -cR01 \
    DB \
    target.fasta

## -v = verbose output
## -S 1 = forward strand only
## -cR01 = lowercase bases in input will be made uppercase, then repetitive regions will be soft-masked with lowercase 

# remove unneeded files
rm -f target.fasta
