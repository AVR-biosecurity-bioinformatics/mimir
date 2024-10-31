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

# params parsing
if ( params.add_root == "true" ) {
    n_columns <- 9
    summary_colnames <- 
        c(
            "Sequences", 
            "tax_ids", 
            "root", 
            "kingdom", 
            "phylum",
            "class", 
            "order", 
            "family", 
            "genus", 
            "species"
        )
} else {
    n_columns <- 8
    summary_colnames <- 
        c(
            "Sequences", 
            "tax_ids", 
            "kingdom", 
            "phylum",
            "class", 
            "order", 
            "family", 
            "genus", 
            "species"
        )
}

### run code

## summarise
taxa_summary <- 
    names(seqs) %>%
    stringr::str_split_fixed(";", n = n_columns) %>%
    tibble::as_tibble() %>%
    tidyr::separate(
        V1, 
        into = c("Sequences", "tax_ids"), 
        sep = "\\|"
    ) %>%
    magrittr::set_colnames(summary_colnames) %>%
    dplyr::summarise_all(n_distinct)

# save as .csv
readr::write_csv(
    x = taxa_summary, 
    file = "taxa_summary.csv"
)

## TODO: generate csv of ids and taxon ranks per sequence
