#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = fasta_file

FASTA_LIST=${3}

# create empty combined output files
touch counts.csv
touch removed.fasta

# loop through all input fasta files, saving retained sequences to lineage-named file
for i in $FASTA_LIST; do

    # get basename of input file
    BASENAME=$( basename $i .fasta )

    # absorb identical sequences and containments, merging names of duplicate sequences into the header of retained sequence (delim: >)
    dedupe.sh \
        mergenames=t \
        absorbrc=f \
        overwrite=t \
        uniquenames=f \
        in=$i \
        out=${BASENAME}.dd.fasta \
        outd=outd.fasta

    # add removed sequences to the bulk removed.fasta file
    cat outd.fasta >> removed.fasta

    # extract just headers from .fasta 
    grep ">" ${BASENAME}.dd.fasta > headers.txt

    # count number of > per line
    sed 's/[^>]//g' headers.txt | awk '{ print length }' > counts.txt

    # get primary sequence header per line
    grep -oE "^>[^>]+" headers.txt > primary.txt

    # paste headers and counts into a csv containing counts for all sequences in the process
    paste primary.txt counts.txt -d , >> counts.csv

    # remove merged names from .fasta and unwrap
    awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ${BASENAME}.dd.fasta | \
        sed 's/^\(>[^>]\+\).*$/\1/' \
        > ${BASENAME}.retained.fasta

done 