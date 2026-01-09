#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = fasta_file

### NOTE: This container (nanozoo/bbmap:38.86--9ebcbfa) uses mawk instead of gawk

FASTA_LIST=${3}

# create empty combined output files
touch counts.tsv
touch removed.fasta

### workflow
# 1. input is a list of genus-level .fasta files
# 2. for each file, use awk to split into N files that each contain a single species lineage
#   a. make sure to save the name of each hash to a temp file that will be used later
#   b. create a genus-specific 'dd.fasta' file that all the retained sequences will be appended to
# 3. for each species file, run BBmap dedupe.sh, with the same code as before, appending removed sequences to a bulk file and appending retained sequences to the genus file

# loop through all input fasta files, saving retained sequences to lineage-named file
for i in $FASTA_LIST; do

    echo "Starting file $i"

    # get basename of input file
    BASENAME=$( basename $i .fasta )

    # make empty file for retained sequences per genus
    touch ${BASENAME}.retained.fasta

    # check input fasta is unwrapped, and if wrapped, unwrap
    if head -n3 $i | grep "^>" | wc -l | grep -q 2; then
        # if first three lines contain two headers, input is unwrapped
        echo "input is unwrapped"
        FASTA=$i
    else 
        # if first three lines don't contain two headers, input is likely wrapped
        # unwrap fasta
        echo "input is wrapped"
        mawk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $i > unwrapped.fasta
        FASTA=unwrapped.fasta
    fi

    # split genus-level file into species-level files, saving sequences to new files based on their lineage string hash
    mawk  -v RS=">" -v ORS="" -v FS="[\r\n]+" -v OFS="\n" ' NR>1 {
        split($1, he, ";" )
        lineage_string=""
        for (i = 2; i <= 8; i++ )
            lineage_string = lineage_string ";" he[i]
        tmp="echo \""lineage_string"\" | md5sum | cut -f1 -d\" \" | tr -d \"\n\""
        tmp | getline linhash
        gsub(/[^[:alnum:]]/, "", linhash)
        print ">" $0 > ( "tmpspp." linhash ".fasta")
        }' \
        < $FASTA

    # loop through temporary species-level files
    for j in tmpspp.*.fasta; do
        # only run dedupe.sh if more than one sequence
        if [[ $( grep -c "^>" $j ) > 1 ]]; then   
            # absorb identical sequences and containments, merging names of duplicate sequences into the header of retained sequence (delim: >)
            dedupe.sh \
                -Xmx6g \
                mergenames=t \
                absorbrc=f \
                overwrite=t \
                uniquenames=f \
                threads=${2} \
                in=$j \
                out=tmp.ret.fasta \
                outd=tmp.rem.fasta

            # add removed sequences to the bulk removed.fasta file
            cat tmp.rem.fasta >> removed.fasta
        else 
            echo "Skipping dedupe.sh for ${j}"
            cp $j tmp.ret.fasta
        fi

        # extract just headers from .fasta 
        grep ">" tmp.ret.fasta > headers.txt

        # count number of > per line
        sed 's/[^>]//g' headers.txt | awk '{ print length }' > counts.txt

        # get primary sequence header per line, without the leading '>'
        grep -oE "^>[^>]+" headers.txt | cut -c2- > primary.txt

        # paste headers and counts into a tsv containing counts for all sequences in the process
        paste primary.txt counts.txt >> counts.tsv

        # remove merged names from .fasta, unwrap and append to genus-level file of retained sequences
        awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' tmp.ret.fasta | \
            sed 's/^\(>[^>]\+\).*$/\1/' \
            >> ${BASENAME}.retained.fasta

        # cleanup tmp files
        rm -f tmp.ret.fasta
        rm -f tmp.rem.fasta
        rm -f headers.txt
        rm -f primary.txt
        rm -f counts.txt

    done

    # remove temp species files
    rm tmpspp.*.fasta

done 
