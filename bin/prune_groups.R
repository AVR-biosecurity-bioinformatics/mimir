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

### run code

#Prune group sizes down to 5, removing all identical sequences first
seqs_pruned <- 
    prune_groups(
        x = seqs, # DNAbin or DNAStringset object
        max_group_size = 5, # max sequences to keep
        discardby = "length", # 'length': discard smallest sequences first; 'random': discard randomly
        dedup = TRUE, # remove sequences with identical taxonomic name and sequence first
        prefer = NULL, # vector of sequence names to prefer (eg. high-quality internal sequences)
        quiet = FALSE
    )

# save filtered sequences as .rds file
saveRDS(seqs_pruned, "seqs_pruned.rds")

# write fasta for debugging
if ( params.all_fasta == "true"){
    write_fasta(
        seqs_pruned, 
        file = "seqs_pruned.fasta", 
        compress = FALSE
        )
}