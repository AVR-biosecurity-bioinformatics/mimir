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
    "taxon"
    )
lapply(nf_vars, nf_var_check)

### process variables 

## TODO: replace with params value
marker <- "COI[GENE] OR COX1[GENE] OR COXI[GENE]"

### run code

## Fetch sequences from GenBank by searching for taxon
if ( stringr::str_detect(taxon, "^\\d+$") ){ # if supplied taxon is a number (uid)...
    genbank_records <- 
        rentrez::entrez_search(
            db = "nuccore",
            term = paste0("(txid", taxon, "[Organism:exp]) AND (",marker,") AND ",params.min_length,":",params.max_length,"[Sequence Length]"),
            retmax = 0
        )
} else {
    genbank_records <- 
        rentrez::entrez_search(
            db = "nuccore",
            term = paste0("(",taxon,"[ORGN]) AND (",marker,") AND ",params.min_length,":",params.max_length,"[Sequence Length]"),
            retmax = 0
        )
}

# get count
records_count <- genbank_records[[2]] 

# write count file
write(x = records_count, file = paste0(records_count %>% format(., scientific = FALSE),".count"))




