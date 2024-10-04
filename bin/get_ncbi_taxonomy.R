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
    "params_dict"
    )
lapply(nf_vars, nf_var_check)

### process variables 



### run code

# download ncbi taxonomy files (can be found in new 'ncbi_taxdump' dir)
ncbi_rankedlineage <- get_ncbi_taxonomy() # import rankedlineage.dmp

# import nodes.dmp
ncbi_nodes <- 
    readr::read_tsv(
        "ncbi_taxdump/nodes.dmp",
        col_names = c(
            "tax_id",
            "parent",
            "rank",
            "embl_code",
            "division_id",
            "div_inherit",
            "gencode",
            "gc_inherit",
            "mtgencode",
            "mtgc_inherit",
            "hidden",
            "hidden_subtree",
            "comments",
            "plgencode",
            "plgc_inherit",
            "specified_species",
            "hdgencode",
            "hdgc_inherit"
        ),
        col_types = c("i-c-c-c-c-i-c-i-c-i-i-i-c-c-i-i-c-i-")
    )

# import fullnamelineage.dmp and keep only taxid and name
ncbi_taxidnames <-
    readr::read_tsv(
        "ncbi_taxdump/fullnamelineage.dmp",
        col_names = c("tax_id", "tax_name", "lineage"),
        col_types = ("i-c-c-")
    ) %>%
    dplyr::select(tax_id, tax_name)

# import names.dmp
ncbi_names <- 
    readr::read_tsv(
        "ncbi_taxdump/names.dmp",
        col_names = c("tax_id","name_txt","unique_name","name_class"),
        col_types = c("i-c-c-c-")
    )

# create tibble with taxid, name and rank
ncbi_taxid_name_rank <- 
    ncbi_taxidnames %>%
    dplyr::left_join(., ncbi_nodes, by = "tax_id") %>%
    dplyr::select(tax_id, tax_name, rank)

# create ncbi synonyms tibble
ncbi_synonyms <- 
    ncbi_names %>%
    dplyr::left_join(., ncbi_taxid_name_rank, by = "tax_id") %>%
    dplyr::filter(name_class == "synonym")

## save files and objects
# save rankedlineage db object
saveRDS(ncbi_rankedlineage, "ncbi_rankedlineage.rds")
# save nodes object
saveRDS(ncbi_nodes, "ncbi_nodes.rds")
# save taxidnames object
saveRDS(ncbi_taxidnames, "ncbi_taxidnames.rds")
# save names object
saveRDS(ncbi_names, "ncbi_names.rds")
# save taxid_name_rank object
saveRDS(ncbi_taxid_name_rank, "ncbi_taxid_name_rank.rds")
# save synonyms object
saveRDS(ncbi_synonyms, "ncbi_synonyms.rds")






