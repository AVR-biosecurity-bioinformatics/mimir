#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = bold_db_path
# $4 = bold_db_url

## check if either parameter is valid
if [[ "$3" = "no_path" && "$4" = "no_url" ]]; then 
    # at least one needs to be specified
    echo "Either 'bold_db_path' or 'bold_db_url' need to be specified!"
    exit 1

# if bold_db_path is a directory
elif [[ -d $3 ]]; then 
    # if uncompressed files exist in directory
    if [ -f ${3}/BOLD_Public*.tsv ] && [ -f ${3}/BOLD_Public.*.datapackage.json ]; then
        # copy uncompressed files to work dir to save as outputs
        cp ${3}/BOLD_Public*.tsv .
        cp ${3}/BOLD_Public.*.datapackage.json . 
    
    # if compressed file exists in directory
    elif [ -f ${3}/BOLD_Public*.tar.gz ]; then
        # decompress compressed database file to current directory
        tar -C . -xzf ${3}/BOLD_Public*.tar.gz 

    else 
        echo "BOLD database files not found in directory $3"
        exit 1
    fi 

# if bold_db_path is a compressed file
elif [[ "$3" =~ .tar.gz$ ]]; then
    # decompress file to current directory
    tar -C . -xzf $3

# if bold_db_url is a URL (hopefully valid)
elif [[ "$4" =~ ^https?:// ]]; then
    # get name of database from URL
    DB_NAME=$( echo "$4" | grep -Po "(?<=id=).+?(?=&uid)" ) # using Perl-style grep so don't have to escape ()
    # download file 
    curl -o "$DB_NAME.tar.gz" "$4"
    # uncompress file
    tar -xzf $DB_NAME.tar.gz

else 
    echo "One or more of bold_db_path ($3) and params.bold_db_url ($4) is invalid."
    exit 1
fi 
