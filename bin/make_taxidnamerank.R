#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    # "Biostrings",
    # "DECIPHER",
    # "IRanges",
    # "R.utils",
    # "RCurl",
    # "ape",
    # "aphid",
    # "bold",
    # "data.table",
    # "data.tree",
    # "dplyr",
    # "entropy",
    # "fs",
    # "furrr",
    # "future",
    # "httr",
    # "kmer",
    # "magrittr",
    # "methods",
    # "openssl",
    # "parallel",
    # "phytools",
    # "purrr",
    # "readr",
    # "rentrez",
    # "rvest",
    # "stats",
    # "stringr",
    # "taxize",
    # "tibble",
    # "tidyr",
    # "utils",
    # "vroom",
    # "xml2",
    "tidyverse",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "params_dict",
    "ncbi_taxidnames",
    "ncbi_nodes"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# import tibbles
ncbi_taxidnames <- readRDS(ncbi_taxidnames)
ncbi_nodes <- readRDS(ncbi_nodes)

### run code

# make tibble of ncbi_taxidnamerank
ncbi_taxidnamerank <- 
    ncbi_taxidnames %>%
    dplyr::left_join(., ncbi_nodes, by = "tax_id") %>%
    dplyr::select(tax_id, tax_name, rank)

# save taxidnamerank db object
saveRDS(ncbi_taxidnamerank, "ncbi_taxidnamerank.rds")

# remove large objects to make saving R environment faster when there are no errors
rm(ncbi_taxidnames)
rm(ncbi_nodes)
rm(ncbi_taxidnamerank)






