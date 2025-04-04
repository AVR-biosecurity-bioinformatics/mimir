#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = fasta_file

# clustalo \
#     -i $3 \
#     -o aligned.fasta \
#     --outfmt=fasta \
#     --wrap=0 \
#     --seqtype=DNA \
#     --dealign \
#     --threads=${2}

### clustalo code
# ## workaround issue of clustalo truncating sequence headers
# # get headers
# awk -v RS=">" -v ORS="\n" -v FS="[\r\n]+" -v OFS="\n" 'NR>1 { print $1 }' $3 > headers.txt
# # truncate headers to seqid + taxid
# cat $3 | sed 's/;.*//' > trunc.fasta
# # align file with truncated headers, keeping original order
# clustalo \
#     -i trunc.fasta \
#     -o aligned.trunc.fasta \
#     --infmt=fasta \
#     --outfmt=fasta \
#     --wrap=999999 \
#     --seqtype=DNA \
#     --dealign \
#     --output-order=input-order \
#     --threads=${2}
# # replace with original headers
# awk 'NR==FNR{a[NR]=$0;next}FNR%2{$0=">"a[int(FNR/2+1)]}1' headers.txt aligned.trunc.fasta > aligned.fasta

### mafft code
# set tmp dir for files
# export MAFFT_TMPDIR="/tmp"

# mafft \
#     --nuc \
#     --thread $2 \
#     --linelength -1 \
#     --auto \
#     "${3}" \
#     > aligned.fasta


### for family-level alignment
export MAFFT_TMPDIR="/tmp"

# if only one sequence in file, skip alignment
if [[ $( grep -c "^>" $3 ) > 1 ]]; then 
    # align
    mafft \
        --nuc \
        --thread $2 \
        --linelength -1 \
        --auto \
        ${3} \
        > aligned.fasta
else 
    # rename as aligned
    cp $3 aligned.fasta
fi 
# add end of line to allow for merging
echo >> aligned.fasta