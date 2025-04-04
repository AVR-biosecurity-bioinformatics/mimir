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
    "db_path"
    )
lapply(nf_vars, nf_var_check)

### process variables 


### run code
# increase timeout 
options(timeout = max(300, getOption("timeout")))

# make tibble of rankedlineage
ncbi_names <- 
    readr::read_tsv(
        paste0(db_path, "/names.dmp"),
        col_names = c("tax_id","name_txt","unique_name","name_class"),
        col_types = c("i-c-c-c-")
    )

# save rankedlineage db object
saveRDS(ncbi_names, "ncbi_names.rds")

# remove large objects to make saving R environment faster when there are no errors
rm(ncbi_names)






