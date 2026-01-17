#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "tidyverse",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "params_dict",
    "ncbi_taxidnamerank",
    "ncbi_names",
    "ncbi_lineageparents"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# import tibbles
ncbi_taxidnamerank <- readRDS(ncbi_taxidnamerank)
ncbi_names <- readRDS(ncbi_names)
ncbi_lineageparents <- readRDS(ncbi_lineageparents)

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

# make tibble of ncbi_synonyms
ncbi_synonyms <- 
    ncbi_names %>%
    dplyr::left_join(., ncbi_taxidnamerank, by = "tax_id") %>%
    dplyr::filter(name_class == "synonym") 

ncbi_filteredsynonyms <- 
    ncbi_synonyms %>% 
    # filter to allowed ranks only
    dplyr::filter(rank %in% allowed_ranks) %>% # only keep synonyms for allowed ranks
    dplyr::left_join( # add parent_taxon to verify synonym 
        ., 
        ncbi_lineageparents %>% dplyr::select(tax_id, parent_taxon, grandparent_taxon),
        by = join_by(tax_id)
    ) %>%
    dplyr::select(synonym = name_txt, tax_name, rank, parent_taxon, grandparent_taxon) # only keep four columns
    ## NOTE: 'synonym' is the taxon name synonymous with the valid taxon name 'tax_name', 'rank' is the rank of the taxon

# save db objects
saveRDS(ncbi_synonyms, "ncbi_synonyms.rds")
saveRDS(ncbi_filteredsynonyms, "ncbi_filteredsynonyms.rds")

# remove large objects to make saving R environment faster when there are no errors
rm(ncbi_taxidnamerank)
rm(ncbi_names)
rm(ncbi_synonyms)
rm(ncbi_filteredsynonyms)
rm(ncbi_lineageparents)






