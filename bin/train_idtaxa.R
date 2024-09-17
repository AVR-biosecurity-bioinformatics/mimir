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
    "seqs_file",
    "db_file"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read sequences from file
seqs <- readRDS(seqs_file)

# read sequences from file
db <- readRDS(db_file)

### run code

# train model
training_set <- 
    train_idtaxa(
        seqs, 
        max_group_size = 10, 
        max_iterations = 3,  
        allow_group_removal = TRUE,  
        get_lineage = TRUE, 
        db = db, 
        quiet = FALSE
    )

# save model
saveRDS(training_set, "idtaxa_model.rds")