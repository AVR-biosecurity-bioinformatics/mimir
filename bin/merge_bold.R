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
    "seq_tibble_list"
    )
lapply(nf_vars, nf_var_check)

### process variables and run code

# do all within a pipe to save memory
merged_tibble <- 
    seq_tibble_list %>%
    stringr::str_extract_all(., pattern = "[^\\s,\\[\\]]+") %>% # convert Groovy list to R list
    unlist() %>%
    lapply(., readRDS) %>% # read in tibbles and store as list 
    dplyr::bind_rows() %>% # merge tibbles
    dplyr::distinct() # remove any duplicate rows

# # read in list of tibbles
# seq_tibble_list <- # convert Groovy to R list format
#     stringr::str_extract_all(seq_tibble_list, pattern = "[^\\s,\\[\\]]+") %>% unlist()

# seq_tibbles <- lapply(seq_tibble_list, readRDS) # read in tibbles and store as list 

# ### run code

# # merge tibbles
# merged_tibble <- 
#     dplyr::bind_rows(seq_tibbles) %>%
#     dplyr::distinct() # make sure there are no duplicated rows

## save subset database as .rds
saveRDS(merged_tibble, "bold_db_merged.rds") ### TODO: make this naming work with multiple taxa as pipeline input
