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
    "seqs_file"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read in seqs from file
seqs <- readRDS(seqs_file)

# get basename
file_basename_noext <- 
    basename(seqs_file) %>%
    stringr::str_remove("_filter_phmm\\.rds$")

### run code

## filter out sequences with stop codons
seqs_filtered <- 
    codon_filter(
        x = seqs, 
        genetic_code = params.genetic_code, 
        tryrc = TRUE, 
        resolve_draws = "majority"
    )

# save filtered sequences as .rds file
if ( !is.null(seqs_filtered)) {
    saveRDS(seqs_filtered, paste0(file_basename_noext,"_filter_stop.rds"))
}

# write fasta for debugging
if ( params.all_fasta == "true" && !is.null(seqs_filtered) ){
    write_fasta(
        seqs_filtered, 
        file = paste0(file_basename_noext,"_filter_stop.fasta"), 
        compress = FALSE
        )
}