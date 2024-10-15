#!/usr/bin/env bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = taxon_id
# $4 = taxon_rank
# $5 = marker
# $6 = params.min_length
# $7 = params.max_length
# $8 = params.use_mito

## process variables

TAXON_ID=${3}

TAXON_RANK=${4}

MARKER=${5}

MIN_LENGTH=${6}

MAX_LENGTH=${7}

USE_MITO=${8}

QUERY_NUC="txid${TAXON_ID}[Organism:exp] AND (${MARKER}) AND ${MIN_LENGTH}:${MAX_LENGTH}[Sequence Length]"

# query by marker name
esearch \
    -db nuccore \
    -query "$QUERY_NUC" \
    | \
    efetch \
    -format acc \
    -mode txt \
    > ${TAXON_ID}_nuc.txt 

# query mitochondrial genomes
if [[ $USE_MITO == "true" ]]; then

    QUERY_MITO="txid${TAXON_ID}[Organism:exp] AND mitochondrion[filter] AND genome AND complete"

    esearch \
        -db nuccore \
        -query "$QUERY_MITO" \
        | \
        efetch \
        -format acc \
        -mode txt \
        > ${TAXON_ID}_mito.txt 

else 
    touch ${TAXON_ID}_mito.txt
fi

# combine accession lists, sort to remove duplicates, then shuffle to mix genomes amongst smaller sequences
cat ${TAXON_ID}_nuc.txt ${TAXON_ID}_mito.txt | \
    sort -u | \
    shuf \
    > ${TAXON_ID}_seqs.txt
