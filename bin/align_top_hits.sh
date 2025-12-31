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

# split input table by query, explicitly sorting by query beforehand
cat $3 | \
    sort -k1,1 - | \
    awk '
        BEGIN {
            prev_query = ""
            file_index = 0
        }
        { 
        if ( $1 != prev_query ) {
            file_index++
        }
        print $0 > file_index".split" 
        prev_query = $1
    }' - 

# for each query, extract the query and all target sequences into a single .fasta file, then align using mafft in a pipe
for i in *.split; do
    # get query
    QUERY=$(head $i -n1 | cut -f1)
    # get query sequence
    grep --no-group-separator -A1 -F "$QUERY" unwrapped.fasta > query.fasta
    # get target sequences
    cat $i | cut -f2 > targets.tmp
    grep --no-group-separator -A1 -F -f targets.tmp unwrapped.fasta > targets.fasta
    # cat into single .fasta, making sure query is first
    cat query.fasta targets.fasta > combined.fasta
    # if list of sequences only contains one, skip alignment
    if [[ $( grep -c ">" combined.fasta ) > 1 ]]; then
        # get sequences that match header, then align with mafft
        cat combined.fasta | \
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
        touch $(basename $i .split).aligned.fasta 
    fi
    # remove temp files
    rm -f *.tmp
    rm -f query.fasta
    rm -f targets.fasta
    rm -f combined.fasta
done

# remove unneeded files
rm -f *.split
rm unwrapped.fasta