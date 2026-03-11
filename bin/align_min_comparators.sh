#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = cfa file

touch alignment.cfa
COUNT=1
while read i; do
    # convert .cfa line to a temp file
    echo "$i" | cat - | sed -e 's/>>>/\n/g' > tmp.fasta
    # align temp file and append to .cfa output file
    mafft \
        --nuc \
        --thread ${2} \
        --linelength -1 \
        --globalpair \
        --maxiterate 1000 \
        --inputorder \
        --quiet \
        tmp.fasta | \
    sed -e ':a;N;$!ba;s/\n/>>>/g' \
    >> alignment.cfa
    echo "Finished alignment $COUNT"
    ((COUNT++))
    rm -f tmp.fasta
done < $3