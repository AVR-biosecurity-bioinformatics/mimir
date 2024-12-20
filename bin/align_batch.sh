#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = fasta_file

# $3 is a list of .fasta files

# for FILE in $3
# do
#     # get basename of file
#     BASE=$( basename "$FILE" .fasta )
#     # do alignment only if number of sequences is >1 AND sequence headers don't contain ";Unclassified" at the end
#     if [[ $( grep -c "^>" $FILE ) > 1 ]] && [[ $( head $FILE -n1 | grep -E ";Unclassified$" | wc -l ) = 0 ]]; then
#         clustalo \
#             -i $FILE \
#             -o ${BASE}.aligned.fasta \
#             --infmt=fasta \
#             --outfmt=fasta \
#             --wrap=999999 \
#             --seqtype=DNA \
#             --dealign \
#             --output-order=input-order \
#             --threads=${2}
#     else 
#         # else rename input file as output file (no alignment needed)
#         cp $FILE ${BASE}.aligned.fasta
#     fi
# done

### workaround clustalo truncating sequence headers >127 characters (which changes lineage string)
## this pre-truncates sequence headers then replaces with the original header after alignment
for FILE in $3
do
    # get basename of file
    BASE=$( basename "$FILE" .fasta )
    # do alignment only if number of sequences is >1 AND sequence headers don't contain ";Unclassified" at the end
    if [[ $( grep -c "^>" $FILE ) > 1 ]] && [[ $( head $FILE -n1 | grep -E ";Unclassified$" | wc -l ) = 0 ]]; then
        ## workaround issue of clustalo truncating sequence headers
        # get headers
        awk -v RS=">" -v ORS="\n" -v FS="[\r\n]+" -v OFS="\n" 'NR>1 { print $1 }' $FILE > ${BASE}.headers.txt
        # truncate headers to seqid + taxid
        cat $FILE | sed 's/;.*//' > ${BASE}.trunc.fasta
        # align file with truncated headers, keeping original order
        clustalo \
            -i ${BASE}.trunc.fasta \
            -o ${BASE}.aligned.trunc.fasta \
            --infmt=fasta \
            --outfmt=fasta \
            --wrap=999999 \
            --seqtype=DNA \
            --dealign \
            --output-order=input-order \
            --threads=${2}
        # replace with original headers
        awk 'NR==FNR{a[NR]=$0;next}FNR%2{$0=">"a[int(FNR/2+1)]}1' ${BASE}.headers.txt ${BASE}.aligned.trunc.fasta > ${BASE}.aligned.fasta
    else 
        # else rename input file as output file (no alignment needed)
        cp $FILE ${BASE}.aligned.fasta
    fi
done



