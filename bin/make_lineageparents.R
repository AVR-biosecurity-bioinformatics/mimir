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
    "ncbi_taxidnamerank"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# import tibbles
ncbi_rankedlineage <- readRDS(ncbi_rankedlineage)
ncbi_taxidnamerank <- readRDS(ncbi_taxidnamerank)

### run code

# make tibble of ncbi_lineageparents
ncbi_lineageparents <- 
    ncbi_rankedlineage %>%
    # remove higher ranks as not used
    dplyr::select(-domain_realm) %>%
    # remove rows where species is not NA (ie. subspecies)
    dplyr::filter(is.na(species)) %>%
    # join to ncbi_taxidnamerank to get rank
    dplyr::left_join(. ,ncbi_taxidnamerank, by = c("tax_id","tax_name")) %>%
    # remove taxa that don't have an allowed rank 
    dplyr::filter(rank %in% c(
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"
        )
    ) %>%
    # mutate to replace the lowest rank with the tax_name (based on 'rank' column value)
    dplyr::mutate(
        species = if_else(rank == "species", tax_name, species),
        genus = if_else(rank == "genus", tax_name, genus),
        family = if_else(rank == "family", tax_name, family),
        order = if_else(rank == "order", tax_name, order),
        class = if_else(rank == "class", tax_name, class),
        phylum = if_else(rank == "phylum", tax_name, phylum),
        kingdom = if_else(rank == "kingdom", tax_name, kingdom)
    ) %>% 
    # add parent_rank, grandparent_rank
    dplyr::mutate(
        parent_rank = case_match(
            rank,
            "species" ~ "genus",
            "genus" ~ "family",
            "family" ~ "order",
            "order" ~ "class",
            "class" ~ "phylum",
            "phylum" ~ "kingdom",
            "kingdom" ~ NA
        ),
        grandparent_rank = case_match(
            parent_rank,
            "species" ~ "genus",
            "genus" ~ "family",
            "family" ~ "order",
            "order" ~ "class",
            "class" ~ "phylum",
            "phylum" ~ "kingdom",
            "kingdom" ~ NA
        ),
        `NA` = NA # need this to allow 'get()' below to work if "*_rank" is NA
    ) %>% 
    dplyr::rowwise() %>%
    # get parent_taxon and grandparent_taxon
    dplyr::mutate(
        parent_taxon = get(parent_rank),
        grandparent_taxon = get(grandparent_rank)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-`NA`)

# save db object
saveRDS(ncbi_lineageparents, "ncbi_lineageparents.rds")

# remove large objects to make saving R environment faster when there are no errors
rm(ncbi_rankedlineage)
rm(ncbi_taxidnamerank)
rm(ncbi_lineageparents)






