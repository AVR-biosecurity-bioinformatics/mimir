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
    "fasta_file"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read in fasta file as DNAbin
seqs <- ape::read.FASTA(fasta_file)


### run code

# get names of sequences (.fasta header)
seq_names <- names(seqs)

# check names are in format "[accession]|[internal taxid];[lineage string]"
name_check <- stringr::str_detect(seq_names, "^[A-Za-z0-9_\\-.]+?\\|[A-Za-z0-9:]+?;\\w+?;\\w+?;\\w+?;\\w+?;\\w+?;\\w+?;[A-Za-z0-9_ .\\-]+?$")

if (!all(name_check)) {
  stop("ERROR: Not all sequence headers in '--internal_seqs' .fasta are formatted correctly.")
}

# replace pure ID with "INTERNAL:ID" (if NCBI: not detected)
seq_names_new <- stringr::str_replace(seq_names, "\\|(?!INTERNAL:|NCBI:)", "\\|INTERNAL:")

# update names
names(seqs) <- seq_names_new

# write .fasta to file
write.FASTA(seqs, file = "internal_imported.fasta")
