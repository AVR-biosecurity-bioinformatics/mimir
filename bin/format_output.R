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
    "fasta_file",
    "add_root",
    "aligned_output",
    "compressed_output"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read fasta from file
seqs <- ape::read.FASTA(fasta_file)

# parse add_root
if ( add_root == "true" ) {
    add_root <- TRUE
} else if ( add_root == "false" ) {
    add_root <- FALSE
} else {
    stop("ERROR: '--add_root' should be 'true' or 'false'.")
}

# parse aligned_output
if ( aligned_output == "true" ) {
    aligned_output <- TRUE
} else if ( aligned_output == "false" ) {
    aligned_output <- FALSE
} else {
    stop("ERROR: '--aligned_output' should be 'true' or 'false'.")
}

# parse compressed_output
if ( compressed_output == "true" ) {
    compressed_output <- TRUE
} else if ( compressed_output == "false" ) {
    compressed_output <- FALSE
} else {
    stop("ERROR: '--compressed_output' should be 'true' or 'false'.")
}

### run code

# add "Root" to taxonomic hierarchy if requested
if ( add_root ) {
    print("Adding 'Root' to taxonomic hierarchy")
    names(seqs) <- 
        names(seqs) %>%
        stringr::str_replace(";",";Root;")
}

if ( !aligned_output ){
    print("Dealigning sequences")
    seqs <- 
        seqs %>% 
        DNAbin2DNAstringset(.) %>% 
        DECIPHER::RemoveGaps(., removeGaps = "all") %>%
        ape::as.DNAbin()
}

### TODO: add options to remove spaces, semicolons etc. in sequence headers

### TODO: put marker and taxon information in the final name of the database

# write fasta 
if (compressed_output) {
    write_fasta(
        seqs, 
        file = "final_database.fasta.gz", 
        compress = TRUE
    )
} else {
    write_fasta(
        seqs, 
        file = "final_database.fasta", 
        compress = FALSE
    )
}
