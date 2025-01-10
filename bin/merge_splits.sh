#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads

## files are staged in their own subdirectory, in the format 'dir*/*', eg. 'dir1/[lineage_string].lineage.fasta'

# from https://unix.stackexchange.com/questions/435794/create-new-concatenated-files-of-same-name-in-multiple-directories

# make new directory for merged files
mkdir merged

# merge all files with the same name
find . -name '*.fasta' \
    -exec sh -c 'for pathname do cat "$pathname" >>"merged/${pathname##*/}"; done' find-sh {} +
