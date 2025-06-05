#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    # "Biostrings",
    # "DECIPHER",
    # "IRanges",
    # "R.utils",
    # "RCurl",
    # "ape",
    # "aphid",
    # "bold",
    # "data.table",
    # "data.tree",
    # "dplyr",
    # "entropy",
    # "fs",
    # "furrr",
    # "future",
    # "httr",
    # "kmer",
    # "magrittr",
    # "methods",
    # "openssl",
    # "parallel",
    # "phytools",
    # "purrr",
    # "readr",
    # "rentrez",
    # "rvest",
    # "stats",
    # "stringr",
    # "taxize",
    # "tibble",
    # "tidyr",
    # "utils",
    # "vroom",
    # "xml2",
    "tidyverse",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "params_dict",
    "ncbi_rankedlineage",
    "ncbi_nodes"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# import tibbles
ncbi_rankedlineage <- readRDS(ncbi_rankedlineage)
ncbi_nodes <- readRDS(ncbi_nodes)

### run code

# make tibble of ncbi_gencodes
ncbi_gencodes <- 
    ncbi_rankedlineage %>%
    # remove higher ranks as not used
    dplyr::select(-tidyselect::any_of(realm, domain)) %>%
    # make lineage information consistent with those found in imported sequences
    dplyr::mutate(
        dplyr::across(species:kingdom, .fns = ~replace(., is.na(.), "Unclassified")), # replace NA with "Unclassified"
        dplyr::across(species:kingdom, .fns = ~stringr::str_replace_all(., " +", " ")) # replace runs of spaces with a single space 
    ) %>%
    # join to nodes data
    dplyr::left_join(., ncbi_nodes, by = "tax_id") %>%
    # retain columns of interest
    dplyr::select(
        tax_id,
        tax_name,
        rank,
        kingdom,
        phylum, 
        class,
        order,
        family,
        genus,
        species,
        gencode, 
        mtgencode,
        plgencode,
        hdgencode
    ) %>%
    # filter to typical ranks
    dplyr::filter(rank %in% c("kingdom","phylum", "class", "order", "family", "genus", "species")) %>%
    # populate lowest rank with name
    dplyr::mutate(
        species = if_else(rank == "species", tax_name, species),
        genus = if_else(rank == "genus", tax_name, genus),
        family = if_else(rank == "family", tax_name, family),
        order = if_else(rank == "order", tax_name, order),
        class = if_else(rank == "class", tax_name, class),
        phylum = if_else(rank == "phylum", tax_name, phylum),
        kingdom = if_else(rank == "kingdom", tax_name, kingdom)
    )

# save rankedlineage db object
saveRDS(ncbi_gencodes, "ncbi_gencodes.rds")

# remove large objects to make saving R environment faster when there are no errors
rm(ncbi_rankedlineage)
rm(ncbi_nodes)
rm(ncbi_gencodes)






