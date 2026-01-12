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
    "alignment_files",
	"seqs_file",
    "rf_counts_tsv",
    "thresholds_csv",
    "con_min_n",
    "con_min_prop"
)
lapply(nf_vars, nf_var_check)

### process variables 

# read in list of sequences
alignment_files <-
    stringr::str_extract_all(alignment_files, pattern = "[^\\s,\\[\\]]+") %>% unlist()

# read in sequences as list
alignment_list <- lapply(alignment_files, Biostrings::readDNAStringSet)

# read in all sequences
seqs <- Biostrings::readDNAStringSet(seqs_file)

# read in rf_counts csv
rf_counts <- readr::read_tsv(rf_counts_tsv, col_names = c("name","n"), show_col_types = FALSE)

# read in thresholds csv
thresholds <- readr::read_csv(thresholds_csv,  show_col_types = FALSE)

# consensus parameters
con_min_n <- as.numeric(con_min_n)
con_min_prop <- as.numeric(con_min_prop)

### run code

stop()


