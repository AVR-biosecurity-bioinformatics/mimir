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

# ranks vector of only ranks we want to keep
allowed_ranks <-
    c(
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"
    )

### run code
# download ncbi taxonomy files (can be found in new 'ncbi_taxdump' dir)
message("Downloading NCBI taxdump")
ncbi_rankedlineage <- get_ncbi_taxonomy() # import rankedlineage.dmp

# import nodes.dmp
message("Importing nodes.dmp")
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
message("Importing fullnamelineage.dmp")
ncbi_taxidnames <-
    readr::read_tsv(
        "ncbi_taxdump/fullnamelineage.dmp",
        col_names = c("tax_id", "tax_name", "lineage"),
        col_types = ("i-c-c-")
    ) %>%
    dplyr::select(tax_id, tax_name)

# import names.dmp
message("Importing names.dmp")
ncbi_names <- 
    readr::read_tsv(
        "ncbi_taxdump/names.dmp",
        col_names = c("tax_id","name_txt","unique_name","name_class"),
        col_types = c("i-c-c-c-")
    )

# create tibble with taxid, name and rank
message("Creating ncbi_taxnamerank")
ncbi_taxidnamerank <- 
    ncbi_taxidnames %>%
    dplyr::left_join(., ncbi_nodes, by = "tax_id") %>%
    dplyr::select(tax_id, tax_name, rank)

# create ncbi synonyms tibble
message("Creating ncbi_synonyms")
ncbi_synonyms <- 
    ncbi_names %>%
    dplyr::left_join(., ncbi_taxidnamerank, by = "tax_id") %>%
    dplyr::filter(name_class == "synonym")

# create 'ncbi_lineageparents' tibble
message("Creating ncbi_lineageparents")
ncbi_lineageparents <- 
    ncbi_rankedlineage %>%
    # remove superkingdom as not used
    dplyr::select(-superkingdom) %>%
    # remove rows where species is not NA (ie. subspecies)
    dplyr::filter(is.na(species)) %>%
    # join to ncbi_taxidnamerank to get rank
    dplyr::left_join(. ,ncbi_taxidnamerank, by = c("tax_id","tax_name")) %>%
    # remove taxa that don't have an allowed rank 
    dplyr::filter(rank %in% allowed_ranks) %>%
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

# create ncbi_rankedlineage_noname
ncbi_rankedlineage_noname <-
    ncbi_rankedlineage %>%
    dplyr::select(-superkingdom) %>%
    dplyr::left_join(., ncbi_taxidnamerank, by = c("tax_id","tax_name")) %>%
    dplyr::mutate(
            species = if_else(rank == "species", tax_name, species),
            genus = if_else(rank == "genus", tax_name, genus),
            family = if_else(rank == "family", tax_name, family),
            order = if_else(rank == "order", tax_name, order),
            class = if_else(rank == "class", tax_name, class),
            phylum = if_else(rank == "phylum", tax_name, phylum),
            kingdom = if_else(rank == "kingdom", tax_name, kingdom)
        ) %>%
    dplyr::select(-tax_name, -rank)

# create ncbi_gencodes (genetic codes per taxid)
ncbi_gencodes <- 
    ncbi_rankedlineage %>%
    # remove superkingdom rank
    dplyr::select(-superkingdom) %>%
    # make lineage information consistent with those found in imported sequences (ie. problematic characters replaced with underscores)
    dplyr::mutate(
        dplyr::across(species:kingdom, .fns = ~replace(., is.na(.), "Unclassified")), # replace NA with "Unclassified"
        dplyr::across(species:kingdom, .fns = ~stringr::str_replace_all(., "[ \\/:\\(\\)&,]", "_")), # replace problematic characters in lineage string with underscores
        dplyr::across(species:kingdom, .fns = ~stringr::str_replace_all(., "_+", "_")) # replace two or more underscores in a row with a single underscore in lineage string
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

## save files and objects
message("Saving .rds files")
# save rankedlineage db object
saveRDS(ncbi_rankedlineage, "ncbi_rankedlineage.rds")
# save nodes object
saveRDS(ncbi_nodes, "ncbi_nodes.rds")
# save taxidnames object
saveRDS(ncbi_taxidnames, "ncbi_taxidnames.rds")
# save names object
saveRDS(ncbi_names, "ncbi_names.rds")
# save taxid_name_rank object
saveRDS(ncbi_taxidnamerank, "ncbi_taxidnamerank.rds")
# save synonyms object
saveRDS(ncbi_synonyms, "ncbi_synonyms.rds")
# save ncbi_lineageparents object
saveRDS(ncbi_lineageparents, "ncbi_lineageparents.rds")
# save ncbi_rankedlineage_noname object
saveRDS(ncbi_rankedlineage_noname, "ncbi_rankedlineage_noname.rds")
# save ncbi_gencodes object
saveRDS(ncbi_gencodes, "ncbi_gencodes.rds")

# remove large objects to make saving R environment faster when there are no errors
rm(ncbi_rankedlineage)
rm(ncbi_nodes)
rm(ncbi_taxidnames)
rm(ncbi_names)
rm(ncbi_taxidnamerank)
rm(ncbi_synonyms)
rm(ncbi_lineageparents)
rm(ncbi_rankedlineage_noname)
rm(ncbi_gencodes)







