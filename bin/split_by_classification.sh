#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = list of fasta files

# create lists of files that are either partially classified or fully classified
touch full.list
touch part.list

for i in ${3}; do
    # if the lineage string ends with two unclassified ranks, the genus is considered partially classified
    if head $i -n1 | grep -q ";Unclassified;Unclassified$" - ; then
        echo $i >> part.list
    # otherwise it is considered fully classified (to genus)
    else
        echo $i >> full.list
    fi
done 

## rename files based on list membership
# partially classified lineages
while IFS= read -r line; do 
    mv "$line" "$(basename $line .fasta).part.fasta" ; 
done < part.list
# fully classified lineages
while IFS= read -r line; do 
    mv "$line" "$(basename $line .fasta).full.fasta" ; 
done < full.list

# output channels group files based on the end of their file name (.part or .full)