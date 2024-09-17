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
    "seqs_list"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read in list of sequences
seqs_list <- # convert Groovy to R list format
    stringr::str_extract_all(seqs_list, pattern = "[^\\s,\\[\\]]+") %>% unlist()

seqs_list <- lapply(seqs_list, readRDS) # read in taxtabs and store as list of tibbles

### run code

## filter out sequences with stop codons
seqs_combined <- 
    concat_DNAbin(
        seqs_list
    )

# # save filtered sequences as .rds file
# saveRDS(seqs_combined, "seqs_combined.rds")

# write sequences to fasta
write_fasta(
    seqs_combined, 
    file = "combined_chunks.fasta", 
    compress = FALSE
    )
