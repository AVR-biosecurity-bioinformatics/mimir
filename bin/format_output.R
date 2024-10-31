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

# read fasta from file
seqs <- ape::read.FASTA(fasta_file)

# parse params.compressed_output
if ( params.compressed_output == "true" ) {
    compressed_output <- TRUE
} else if ( params.compressed_output == "false" ) {
    compressed_output <- FALSE
} else {
    stop("ERROR: '--compressed_output' should be 'true' or 'false'.")
}

### run code

# add "Root" to taxonomic hierarchy if requested
if ( params.add_root == "true" ) {
    names(seqs) <- 
        names(seqs) %>%
        stringr::str_replace(";",";Root;")
}

### TODO: add options to remove spaces, semicolons etc. in sequence headers

### TODO: put marker and taxon information in the final name of the database

# save renamed sequences as .rds file
saveRDS(seqs, "final_database.rds")

# write fasta 
write_fasta(
    seqs, 
    file = "final_database.fasta", 
    compress = compressed_output
    )