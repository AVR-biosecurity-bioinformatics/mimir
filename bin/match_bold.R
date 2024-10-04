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
    # "merged_tibble",
    "seq_tibble_list",
    "ncbi_rankedlineage",
    "ncbi_taxidnamerank",
    "ncbi_synonyms"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# merged_tibble <- readRDS(merged_tibble)

# create bold_db_targets (do in pipe to save memory)
bold_db_targets <- 
    seq_tibble_list %>%
    stringr::str_extract_all(., pattern = "[^\\s,\\[\\]]+") %>% # convert Groovy list to R list
    unlist() %>%
    lapply(., readRDS) %>% # read in tibbles and store as list 
    dplyr::bind_rows() %>% # merge tibbles
    dplyr::distinct() %>% # remove any duplicate rows
    dplyr::select( # select only needed columns
        processid,
        bold_taxid,
        kingdom,
        phylum,
        class,
        order,
        family,
        genus,
        species,
        nuc
    )            

# read in ncbi db files
ncbi_rankedlineage <-   readRDS(ncbi_rankedlineage)
ncbi_taxidnamerank <-   readRDS(ncbi_taxidnamerank)
ncbi_synonyms <-        readRDS(ncbi_synonyms)

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

## modifying ncbi objects
ncbi_taxdata <- 
    ncbi_rankedlineage %>%
    # remove superkingdom as not used
    dplyr::select(-superkingdom) %>%
    # remove rows where species is not NA (which means taxon is subspecies)
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

# create a filtered synonyms tibble
ncbi_synonyms_filtered <- 
    ncbi_synonyms %>% 
    dplyr::filter(rank %in% allowed_ranks) %>% # only keep synonyms for allowed ranks
    dplyr::select(synonym = name_txt, tax_name, rank) # only keep three columns
    ## NOTE: 'synonym' is the taxon name synonymous with the valid taxon name 'tax_name', 'rank' is the rank of the taxon

### TODO: manually add in known synonyms between BOLD and NCBI
## for example: family "Xylophagaidae" in NCBI is called "Xylophagidae" in BOLD


## replace BOLD taxon name with NCBI name if they're the same rank and an NCBI synonym
# do this while tracking replacement - create columns that contain old names if they're replaced
bold_db_targets_syns <- 
    bold_db_targets %>%
    dplyr::mutate(across(kingdom:species, .fns = ~replace(., . == "None", NA))) %>% # replace "None" with NA within rank columns
    ## for each rank, find and replace based on ncbi synonyms
    # kingdom
    dplyr::left_join(., ncbi_synonyms_filtered %>% dplyr::filter(rank == "kingdom"), by = join_by(kingdom == synonym)) %>%
    dplyr::mutate(kingdom_old = if_else(!is.na(tax_name), kingdom, NA), kingdom = if_else(!is.na(tax_name), tax_name, kingdom)) %>% 
    dplyr::select(-tax_name, -rank) %>%
    # phylum
    dplyr::left_join(., ncbi_synonyms_filtered %>% dplyr::filter(rank == "phylum"), by = join_by(phylum == synonym)) %>%
    dplyr::mutate(phylum_old = if_else(!is.na(tax_name), phylum, NA), phylum = if_else(!is.na(tax_name), tax_name, phylum)) %>% 
    dplyr::select(-tax_name, -rank) %>%
    # class
    dplyr::left_join(., ncbi_synonyms_filtered %>% dplyr::filter(rank == "class"), by = join_by(class == synonym)) %>%
    dplyr::mutate(class_old = if_else(!is.na(tax_name), class, NA), class = if_else(!is.na(tax_name), tax_name, class)) %>% 
    dplyr::select(-tax_name, -rank) %>%
    # order
    dplyr::left_join(., ncbi_synonyms_filtered %>% dplyr::filter(rank == "order"), by = join_by(order == synonym)) %>%
    dplyr::mutate(order_old = if_else(!is.na(tax_name), order, NA), order = if_else(!is.na(tax_name), tax_name, order)) %>% 
    dplyr::select(-tax_name, -rank) %>%
    # family
    dplyr::left_join(., ncbi_synonyms_filtered %>% dplyr::filter(rank == "family"), by = join_by(family == synonym)) %>%
    dplyr::mutate(family_old = if_else(!is.na(tax_name), family, NA), family = if_else(!is.na(tax_name), tax_name, family)) %>% 
    dplyr::select(-tax_name, -rank) %>%
    # genus
    dplyr::left_join(., ncbi_synonyms_filtered %>% dplyr::filter(rank == "genus"), by = join_by(genus == synonym)) %>%
    dplyr::mutate(genus_old = if_else(!is.na(tax_name), genus, NA), genus = if_else(!is.na(tax_name), tax_name, genus)) %>% 
    dplyr::select(-tax_name, -rank) %>%
    # species
    dplyr::left_join(., ncbi_synonyms_filtered %>% dplyr::filter(rank == "species"), by = join_by(species == synonym)) %>%
    dplyr::mutate(species_old = if_else(!is.na(tax_name), species, NA), species = if_else(!is.na(tax_name), tax_name, species)) %>% 
    dplyr::select(-tax_name, -rank) %>%
    ## if species was replaced but not the genus, use species name to correct genus name
    dplyr::mutate(
        genus_old = if_else(
            !stringr::str_starts(species, genus) & !is.na(species), # if genus and species don't align and species isn't NA...
            genus, # record old genus name
            NA # else keep NA
        ),
        genus = if_else(
            is.na(genus_old), # if genus_old wasn't populated above...
            genus, # keep original genus name
            stringr::str_extract(species, "^\\w+(?= )") # else use start of species binomial as new genus name
        )
    )
    ### TODO: Do some proper testing of the above code to check it is handling the resolution correctly - check the "synchanges.csv" output file


# record which taxa names were changed from BOLD to NCBI
### TODO: add explicit counting of these as a summary table
bold_ncbi_synchanges <- 
    bold_db_targets_syns %>%
    dplyr::select(kingdom:species, kingdom_old:species_old) %>%
    dplyr::distinct()

# remove columns that record changes
bold_db_targets_nosyn <- 
    bold_db_targets_syns %>%
    dplyr::select(-dplyr::ends_with("_old"))

### matching BOLD taxonomy to NCBI taxonomy

## create new 'identification_rank' (lowest ID rank) and 'identification' (lowest ID taxa) using synonym-adjusted tibble
bold_id_tibble <- 
    bold_db_targets_nosyn %>%
    dplyr::mutate(
        # create new identification_rank
        identification_rank = purrr::pmap_chr( # map collection of vectors, outputting character class
            .l = dplyr::select(., all_of(allowed_ranks)), # select rank columns
            # function to apply
            .f = ~{
                tax <- c(...) #vector of all input values
                names(tax) <- allowed_ranks # make names of vector the ranks
                revtax <- rev(tax) # reverse the vector (so starts with 'species')
                keepvec <- which.max(!is.na(revtax)) # 
                return(names(keepvec))  # Return the rank name
            }
        ),
        # get the rank above the id rank (parent rank)
        parent_rank = purrr::map(
            .x = identification_rank, 
            .f = ~{
                ind <- match(.x, allowed_ranks)
                return(allowed_ranks[ind - 1])
            }
        ) %>% unlist,
        # get grandparent rank
        grandparent_rank = purrr::map(
            .x = identification_rank, 
            .f = ~{
                ind <- match(.x, allowed_ranks)
                return(allowed_ranks[ind - 2])
            }
        ) %>% unlist
    ) %>%
    dplyr::rowwise() %>% # needed for following 'get()' mutates
    # get the taxon names from the columns named in the "*_rank" columns
    dplyr::mutate( 
        identification = get(identification_rank),
        parent = get(parent_rank),
        grandparent = get(grandparent_rank)
    ) %>%  
    dplyr::ungroup()

# join to ncbi_taxdata based on the three identification taxa and their ranks (replace BOLD with NCBI lineage)
bold_ncbi_joined <- 
    bold_id_tibble %>%
    dplyr::left_join(
        ., 
        ncbi_taxdata, 
        by = join_by(
            identification_rank == rank, 
            identification == tax_name,
            parent == parent_taxon,
            grandparent == grandparent_taxon,
            parent_rank,
            grandparent_rank
        )
    ) %>%
    dplyr::select(!dplyr::ends_with(".x")) %>% # remove taxa columns from the original BOLD data
    dplyr::rename_with(., ~str_remove(., '\\.y$'), .cols = dplyr::ends_with(".y")) # rename taxa columns from NCBI
  
# get only sequences that matched to NCBI taxa
bold_ncbi_matched <- 
    bold_ncbi_joined %>%
    dplyr::filter(!is.na(tax_id))

# get tibble of boldids and the corresponding NCBI taxid, with rank
matching_taxids <- 
    bold_ncbi_matched %>%
    dplyr::select(
        taxon_name = identification, 
        bold_taxid, 
        ncbi_taxid = tax_id, 
        rank = identification_rank
    ) %>%
    dplyr::distinct() # keep only one row per taxon

# get processids of sequences that didn't match the NCBI taxonmy
unmatched_processids <- 
    bold_ncbi_joined %>%
    dplyr::filter(is.na(tax_id)) %>%
    dplyr::pull(processid)

# get unmatched sequences with synonym-resolved BOLD taxonomy
bold_ncbi_unmatched <- 
    bold_id_tibble %>%
    dplyr::filter(processid %in% unmatched_processids)

## give matched and unmatched tibbles a consistent format
# format unmatched
bold_unmatched_seqs <-
    bold_ncbi_unmatched %>% 
    dplyr::select(
        seqid = processid,
        taxid = bold_taxid,
        kingdom:species,
        nuc
    ) %>%
    dplyr::mutate(
        taxid = stringr::str_replace(taxid, "^", "BOLD:"),
        across(kingdom:species, .fns = ~replace(., is.na(.), "Unclassified"))
    )

# format matched
bold_matched_seqs <- 
    bold_ncbi_matched %>% 
    dplyr::select(
        seqid = processid,
        taxid = tax_id,
        kingdom:species,
        nuc
    ) %>%
    dplyr::mutate(
        taxid = stringr::str_replace(taxid, "^", "NCBI:"),
        across(kingdom:species, .fns = ~replace(., is.na(.), "Unclassified"))
    )

# combine matched and unmatched
bold_seqs <- 
  dplyr::bind_rows(bold_matched_seqs, bold_unmatched_seqs)

## create FASTA format from sequence tibble
bold_seqs_prefasta <- 
  bold_seqs %>%
  tidyr::unite("id", c(seqid, taxid), sep = "|") %>% # combine ids into a single column
  tidyr::unite("ranks", kingdom:species, sep = ";") %>% # combine ranks into a single column
  tidyr::unite("header", id:ranks, sep = ";") %>% # create header column
  dplyr::mutate(header = stringr::str_replace(header, "^", ">")) # add ">" to start of header
  
bold_seqs_header <- bold_seqs_prefasta %>% dplyr::pull(header) # get header as vector
bold_seqs_nuc <- bold_seqs_prefasta %>% dplyr::pull(nuc) # get sequence as vector

# interleave into .fasta format 
bold_seqs_fasta <- rbind(bold_seqs_header, bold_seqs_nuc) # create matrix
attributes(bold_seqs_fasta) <- NULL # unset dim 

### save outputs
# tibble of changes to BOLD sequence taxonomy by matching to NCBI synonyms
readr::write_csv(bold_ncbi_synchanges, "synchanges.csv")
# tibble of taxa (from sequences) that can be found in BOLD and NCBI 
readr::write_csv(matching_taxids, "matching_taxids.csv")
# write FASTA of sequences in the correct format
readr::write_lines(bold_seqs_fasta, "bold_seqs.fasta")

# remove merged_tibble object as it can be large and will slow down R environment saving in debug mode
rm(merged_tibble)