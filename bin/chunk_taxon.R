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
    "entrez_key"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# set API key if provided
if ( entrez_key != "no_key" ){
    Sys.setenv(ENTREZ_KEY = entrez_key)
}

## convert taxon name to uid, if necessary
if ( stringr::str_detect(taxon, "^\\d+$") ) { # if already a number, leave as is
    # ensure numeric class
    id <- as.numeric(taxon)
    # get character taxon name from id
    taxon_name <- taxize::id2name(id, db = "ncbi")[[1]] %>% pull(name)
} else {
    # convert taxon name to uid
    id <- 
        taxize::get_uid(
            sci_com = taxon,
            modifier = "Scientific Name"
        )[1]
    # save taxon name variable
    taxon_name <- taxon
}

if ( is.na(id) ){
    stop( paste0("ERROR: Target taxon '",taxon,"' could not be unambiguously converted to UID/taxid.\n\tConsider supplying UID directly."))
}

### run code

## check if target taxon rank is higher than the chunk rank
# get classification of target rank
taxon_classification <- 
    taxize::classification(
        sci_id = id,
        db = "ncbi"
    )
# get actual rank and upstream ranks of target taxon
taxon_ranks <- taxon_classification[[1]] %>% pull(rank)

## if chunk rank is lower than taxon rank, chunk taxon
if ( !is.element(params.chunk_rank, taxon_ranks ) ){
    ## Fetch downstream taxa at a given rank
    taxon_chunked <- 
        taxize::ncbi_downstream(
            id = id, 
            downto = params.chunk_rank
        )

    # subset to just a vector of names
    chunk_vector <- taxon_chunked$childtaxa_id    
} else { 
    ## if chunk rank same or higher than taxon rank, just use taxon as the output
    chunk_vector <- taxon_name
} 

# check that a vector was produced
if ( is.null(chunk_vector) ){
    stop (paste0("ERROR: No taxon list produced from taxid = '",id,"' and chunk_rank = '",params.chunk_rank,"'.\n\tCheck target taxon is not at or below the chunk rank."))
}

# save vector as text file
write(chunk_vector, "tax_list.txt")
