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

# read ncbi db file
db <- readRDS(db_file)

### run code

# reformat to complete taxonomic heirarchy 
seqs_renamed <- 
    taxreturn::reformat_hierarchy(
        x = seqs, 
        db = db, 
        ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
        quiet = FALSE
    )

# add "Root" to the hierarchy
names(seqs_renamed) <- 
  names(seqs_renamed) %>%
  stringr::str_replace(";",";Root;")

# save renamed sequences as .rds file
saveRDS(seqs_renamed, "seqs_renamed.rds")

# write fasta for debugging
if ( params.all_fasta == "true"){
    write_fasta(
        seqs_renamed, 
        file = "seqs_renamed.fasta", 
        compress = FALSE
        )
}