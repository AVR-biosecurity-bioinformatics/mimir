#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "Biostrings",
    "DECIPHER",
    "IRanges",
    "R.utils",
    "RCurl",
    "ape",
    "aphid",
    "bold",
    "data.table",
    "data.tree",
    "dplyr",
    "entropy",
    "fs",
    "furrr",
    "future",
    "httr",
    "kmer",
    "magrittr",
    "methods",
    "openssl",
    "parallel",
    "phytools",
    "purrr",
    "readr",
    "rentrez",
    "rvest",
    "stats",
    "stringr",
    "taxize",
    "tibble",
    "tidyr",
    "utils",
    "vroom",
    "xml2",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "params_dict",
    "taxon",
    "entrez_key",
    "chunk_rank"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# set API key if provided
if ( entrez_key != "no_key" ){
    Sys.setenv(ENTREZ_KEY = entrez_key)
}

# convert taxon to uid, if 
if ( stringr::str_detect(taxon, "^\\d+$") ) {
    id <- as.numeric(taxon)
} else {
    id <- 
        taxize::get_uid(
            sci_com = taxon,
            modifier = "Scientific Name"
        )[1]
}

if ( is.na(id) ){
    stop( paste0("ERROR: Target taxon '",id,"' could not be unambiguously converted to UID/taxid.\n\tConsider supplying UID directly."))
}

### run code

## Fetch downstream taxa at a given rank
taxon_chunked <- 
    taxize::ncbi_downstream(
        id = id, 
        downto = chunk_rank
    )

# subset to just a vector of names
chunk_list <- taxon_chunked$childtaxa_name

# check that a list was produced
if ( is.null(chunk_list) ){
    stop (paste0("ERROR: No taxon list produced from taxid = '",id,"' and chunk_rank = '",chunk_rank,"'.\n\tCheck target taxon is not at or below the chunk rank."))
}

# save list as text file
write(chunk_list, "tax_list.txt")
