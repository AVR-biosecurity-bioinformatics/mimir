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
    "ncbi_taxidnamerank",
    "ncbi_names"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# import tibbles
ncbi_taxidnamerank <- readRDS(ncbi_taxidnamerank)
ncbi_names <- readRDS(ncbi_names)

### run code

# make tibble of ncbi_synonyms
ncbi_synonyms <- 
    ncbi_names %>%
    dplyr::left_join(., ncbi_taxidnamerank, by = "tax_id") %>%
    dplyr::filter(name_class == "synonym")

# save db object
saveRDS(ncbi_synonyms, "ncbi_synonyms.rds")

# remove large objects to make saving R environment faster when there are no errors
rm(ncbi_taxidnamerank)
rm(ncbi_names)
rm(ncbi_synonyms)






