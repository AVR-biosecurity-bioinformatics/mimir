#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = original file name

## files are staged in their own subdirectory, in the format 'dir*/*', eg. 'dir1/[lineage_string].lineage.fasta'

# # from https://unix.stackexchange.com/questions/435794/create-new-concatenated-files-of-same-name-in-multiple-directories

# make new directory for merged files and delete it if it exists already
rm -rf merged
mkdir -p merged

# merge all files with the same name
find . -name '*.fasta' -and -path './dir*' \
    -exec sh -c 'for pathname do cat "$pathname" >>"merged/${pathname##*/}"; done' find-sh {} +


# mkdir -p merged 
# cat dir*/*.fasta > merged/${3}
