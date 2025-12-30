#!/bin/bash
set -e
set -u
set -o pipefail
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = top hits tsv
# $4 = seqs fasta

# unwrap fasta file
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $4 > unwrapped.fasta

# split input table by query, explicitly sorting by query
cat $3 | \
    sort -k1,1 - | \
    awk '
        BEGIN {
            prev_lin = ""
            file_index = 0
        }
        { 
        # get genus lineage string of query
        split($1, query_array, ";")
        lineage_string=""
        for (i = 2; i <= 7; i++ )
            lineage_string = lineage_string ";" query_array[i]
        if (lineage_string != prev_lin) {
            file_index++
        }
        print $0 > file_index".split" 
        prev_lin = lineage_string
    }' - 

# for each query, extract the query and all target sequences into a single .fasta file, then align using mafft in a pipe
for i in *.split; do
    # get query
    QUERY=$(head $i -n1 | cut -f1)
    # get targets
    cat $i | cut -f2 > targets.tmp
    # cat into list
    echo $QUERY | cat - targets.tmp > list.tmp
    # if list of sequences only contains one, skip alignment
    if [[ $( cat list.tmp | wc -l ) > 1 ]]; then
        # get sequences that match header, then align with mafft
        grep --no-group-separator -A1 -F -f list.tmp unwrapped.fasta | \
            mafft \
                --nuc \
                --linelength -1 \
                --globalpair \
                --maxiterate 1000 \
                --thread ${2} \
                --quiet \
                - \
            > $(basename $i .split).aligned.fasta 
        # throw error if output file is empty
        if [ -s $(basename $i .split).aligned.fasta ]; then
            echo "Finished aligning $i"
        else 
            echo "Alignment output for $i is empty"
            exit 1
        fi
    else 
        echo "Skipping alignment for $i"
        touch $(basename $i .split).aligned.fasta 
    fi
    # remove temp files
    rm -f *.tmp
done

# remove unneeded files
rm -f *.split
rm unwrapped.fasta