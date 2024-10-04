#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = gene
# $4 = taxon


datasets summary gene \
    symbol $3 \
    --taxon $4 \
    --as-json-lines \
    | \
    dataformat tsv gene \
    --fields symbol,gene-id,synonyms




