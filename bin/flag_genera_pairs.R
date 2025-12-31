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
    "fasta_files",
    "thresholds_file",
    "seqs_file",
    "counts_file"
)
lapply(nf_vars, nf_var_check)

### process variables 

# read in list of sequences
seqs_list <- 
    stringr::str_extract_all(fasta_files, pattern = "[^\\s,\\[\\]]+") %>% 
    unlist() %>%
    lapply(., Biostrings::readDNAStringSet)

thresholds <- readr::read_csv(thresholds_file, show_col_types = F)

seqs <- Biostrings::readDNAStringSet(seqs_file)

counts <- readr::read_tsv(counts_file, col_names = c("name","n"), show_col_types = FALSE)

stop()


### run code

