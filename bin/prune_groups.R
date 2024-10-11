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

# read sequences from file
seqs <- readRDS(seqs_file)

## parse params
max_group_size <- as.numeric(params.max_group_size)

prune_method <- params.prune_method

### run code

## prune groups of sequences with identical taxonomic IDs down to a certain number
seqs_pruned <- 
    prune_groups(
        x = seqs, # DNAbin or DNAStringset object
        max_group_size = max_group_size, # max sequences to keep
        discardby = prune_method, # 'length': discard smallest sequences first; 'random': discard randomly
        dedup = TRUE, # remove sequences with identical taxonomic name and sequence first
        prefer = NULL, # vector of sequence names to prefer (eg. high-quality internal sequences)
        quiet = FALSE
    )

# save filtered sequences as .rds file
saveRDS(seqs_pruned, "seqs_pruned.rds")

# write fasta
write_fasta(
    seqs_pruned, 
    file = "seqs_pruned.fasta", 
    compress = FALSE
)
