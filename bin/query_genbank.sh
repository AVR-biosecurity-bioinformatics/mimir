#!/usr/bin/env bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = taxon_name
# $4 = taxon_rank
# $5 = marker
# $6 = params.min_length
# $7 = params.max_length

## process variables

TAXON_NAME=${3}

TAXON_RANK=${4}

MARKER=${5}

MIN_LENGTH=${6}

MAX_LENGTH=${7}

QUERY="${TAXON_NAME}[ORGN] AND (${MARKER}) AND ${MIN_LENGTH}:${MAX_LENGTH}[Sequence Length]"

esearch \
    -db nuccore \
    -query "$QUERY" \
    | \
    efetch \
    -format acc \
    -mode txt \
    > ${TAXON_NAME}.txt 