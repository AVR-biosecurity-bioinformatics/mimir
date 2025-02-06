#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = fasta_file
# $4 = primers


##

# check only one sequence in primer file
if [[ $( grep -c "^>" $4 ) > 1 ]]; then
    echo "More than one primer sequence in primer file"
    exit 1
fi

# get name of primer sequence
PRIMER_NAME=$( sed -n 's/^>//p' $4 )

# align primers to database alignment, then retain only primer sequence in output (to save space)
mafft \
    --nuc \
    --thread $2 \
    --linelength -1 \
    --6merpair \
    --keeplength \
    --addfragments $4 \
    "${3}" | \
    awk -vRS=">" -vORS="" -vFS="\n" -vOFS="\n" '
    NR>1 && $1 ~ /^'$PRIMER_NAME'$/ {print ">"$0}
  ' - > aligned_primer.fasta

# throw error if output file is empty
if [ -s aligned_primer.fasta ]; then
    echo "Finished aligning primer sequence"        
else 
    echo "alignment output file is empty"
    exit 1
fi