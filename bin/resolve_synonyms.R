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
    "db_path"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read sequences from file
seqs <- readRDS(seqs_file)

# get basename
file_basename_noext <- 
    basename(seqs_file) %>%
    stringr::str_remove("_filter_(stop|phmm)\\.rds$")

### run code

## filter out sequences with stop codons
seqs_resolved <- 
    resolve_synonyms_ncbi(
        x = seqs,
        dir = db_path
    )

# save filtered sequences as .rds file
saveRDS(seqs_resolved, paste0(file_basename_noext,"_seqs_resolved.rds"))

# write fasta for debugging
if ( params.all_fasta == "true"){
    write_fasta(
        seqs_resolved, 
        file = paste0(file_basename_noext, "_seqs_resolved.fasta"), 
        compress = FALSE
        )
}