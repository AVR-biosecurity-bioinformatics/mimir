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
    "seq_tibble",
    "ncbi_lineageparents",
    "ncbi_filteredsynonyms",
    "genbank_accessions",
    "placeholder_as_unclassified",
    "digits_as_unclassified",
    "bold_idmethod_filter"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# create bold_db_targets (do in pipe to save memory)
bold_db_targets <- 
    readr::read_tsv(seq_tibble) %>%
    dplyr::select(
        processid,
        bold_taxid = taxid,
        kingdom,
        phylum,
        class,
        order,
        family,
        genus, 
        species,
        nuc,
        identification_method,
        insdc_acs, 
        sequence_run_site
    )

# read in ncbi db files
ncbi_lineageparents <- readRDS(ncbi_lineageparents)
ncbi_filteredsynonyms <- readRDS(ncbi_filteredsynonyms)

# read in accessions as vector
gb_acc <- readr::read_lines(genbank_accessions)

# parse if genbank was being used
if (all(gb_acc == c("NOFILE"))){
	use_genbank <- FALSE
} else {
	use_genbank <- TRUE
}

# ranks vector of only ranks we want to keep
allowed_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

# parse params.placeholder_as_unclassified
if ( placeholder_as_unclassified == "true" ){
    placeholder_as_unclassified <- TRUE
} else if ( placeholder_as_unclassified == "false" ){
    placeholder_as_unclassified <- FALSE
} else {
    stop("'placeholder_as_unclassified' must be 'true' or 'false'")
}

# parse params.digits_as_unclassified
if ( digits_as_unclassified == "true" ){
    digits_as_unclassified <- TRUE
} else if ( digits_as_unclassified == "false" ){
    digits_as_unclassified <- FALSE
} else {
    stop("'digits_as_unclassified' must be 'true' or 'false'")
}

# parse params.bold_idmethod_filter
if ( bold_idmethod_filter == "true" ){
    bold_idmethod_filter <- TRUE
} else if ( bold_idmethod_filter == "false" ){
    bold_idmethod_filter <- FALSE
} else {
    stop("'bold_idmethod_filter' must be 'true' or 'false'")
}

### run code

## replace BOLD taxon name with NCBI name if they're the same rank, an NCBI synonym and parent_taxon name also matches (should still work if both are synonyms as gets replaced progressively down the hierarchy)
# do this while tracking replacement - create columns that contain old names if they're replaced
bold_db_targets_syns <- 
    bold_db_targets %>%
    dplyr::mutate(across(kingdom:species, .fns = ~replace(., . == "None", NA))) %>% # replace "None" with NA within rank columns
    ## for each rank, find and replace based on ncbi synonyms
    # kingdom
    dplyr::left_join(., ncbi_filteredsynonyms %>% dplyr::filter(rank == "kingdom"), by = join_by(kingdom == synonym)) %>%
    dplyr::mutate(kingdom_old = if_else(!is.na(tax_name), kingdom, NA), kingdom = if_else(!is.na(tax_name), tax_name, kingdom)) %>% 
    dplyr::select(-tax_name, -rank, -parent_taxon, -grandparent_taxon) %>%
    # phylum
    dplyr::left_join(., ncbi_filteredsynonyms %>% dplyr::filter(rank == "phylum"), by = join_by(phylum == synonym, kingdom == parent_taxon)) %>%
    dplyr::mutate(phylum_old = if_else(!is.na(tax_name), phylum, NA), phylum = if_else(!is.na(tax_name), tax_name, phylum)) %>% 
    dplyr::select(-tax_name, -rank, -grandparent_taxon) %>%
    # class
    dplyr::left_join(., ncbi_filteredsynonyms %>% dplyr::filter(rank == "class"), by = join_by(class == synonym, phylum == parent_taxon, kingdom == grandparent_taxon)) %>%
    dplyr::mutate(class_old = if_else(!is.na(tax_name), class, NA), class = if_else(!is.na(tax_name), tax_name, class)) %>% 
    dplyr::select(-tax_name, -rank) %>%
    # order
    dplyr::left_join(., ncbi_filteredsynonyms %>% dplyr::filter(rank == "order"), by = join_by(order == synonym, class == parent_taxon, phylum == grandparent_taxon)) %>%
    dplyr::mutate(order_old = if_else(!is.na(tax_name), order, NA), order = if_else(!is.na(tax_name), tax_name, order)) %>% 
    dplyr::select(-tax_name, -rank) %>%
    # family
    dplyr::left_join(., ncbi_filteredsynonyms %>% dplyr::filter(rank == "family"), by = join_by(family == synonym, order == parent_taxon, class == grandparent_taxon)) %>%
    dplyr::mutate(family_old = if_else(!is.na(tax_name), family, NA), family = if_else(!is.na(tax_name), tax_name, family)) %>% 
    dplyr::select(-tax_name, -rank) %>%
    # genus
    dplyr::left_join(., ncbi_filteredsynonyms %>% dplyr::filter(rank == "genus"), by = join_by(genus == synonym, family == parent_taxon, order == grandparent_taxon)) %>%
    dplyr::mutate(genus_old = if_else(!is.na(tax_name), genus, NA), genus = if_else(!is.na(tax_name), tax_name, genus)) %>% 
    dplyr::select(-tax_name, -rank) %>%
    # species
    dplyr::left_join(., ncbi_filteredsynonyms %>% dplyr::filter(rank == "species"), by = join_by(species == synonym, genus == parent_taxon, family == grandparent_taxon)) %>%
    dplyr::mutate(species_old = if_else(!is.na(tax_name), species, NA), species = if_else(!is.na(tax_name), tax_name, species)) %>% 
    dplyr::select(-tax_name, -rank) %>%
    dplyr::mutate(across(kingdom:species, .fns = ~replace(., is.na(.), "UNCLASSIFIED"))) %>% # replace NA with UNCLASSIFIED within rank columns (temporarily)
    ## if species was replaced but not the genus, use species name to correct genus name
    dplyr::mutate(
        genus_old = if_else(
            species != "UNCLASSIFIED" & !stringr::str_starts(species, genus), # if genus and species don't align and species isn't NA...
            genus, # record old genus name
            NA # else keep NA
        ),
        genus = if_else(
            is.na(genus_old), # if genus_old wasn't populated above...
            genus, # keep original genus name
            stringr::str_extract(species, "^\\w+(?= )") # else use start of species binomial as new genus name
        )
    ) %>%
    dplyr::mutate(across(kingdom:species, .fns = ~replace(., . == "UNCLASSIFIED", NA))) # replace UNCLASSIFIED with NA within rank columns

# check no rows have been added or removed
if (nrow(bold_db_targets_syns) != nrow(bold_db_targets)) {
    stop(paste0("ERROR: Number of rows in 'bold_db_targets' (",nrow(bold_db_targets),") and 'bold_db_targets_syns (",nrow(bold_db_targets_syns),") doesn't match."))
}

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
                names(tax) <- allowed_ranks # make the names of vector the ranks
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
                return(allowed_ranks[max(ind - 1,1)])
            }
        ) %>% unlist(),
        # get grandparent rank
        grandparent_rank = purrr::map(
            .x = identification_rank,
            .f = ~{
                ind <- match(.x, allowed_ranks)
                return(allowed_ranks[max(ind - 2,1)])
            }
        ) %>% unlist()
    ) %>% 
    # handle cases where identification_rank, parent_rank and/or grand_parent rank are the same
    # this happens in situations where the identification_rank is phylum or kingdom
    dplyr::mutate(
        grandparent_rank = dplyr::if_else(grandparent_rank == parent_rank, NA, grandparent_rank),
        parent_rank = dplyr::if_else(parent_rank == identification_rank, NA, parent_rank)
    ) %>% 
    dplyr::rowwise() %>% # needed for following 'get()' mutates
    # get the taxon names from the columns named in the "*_rank" columns
    dplyr::mutate( 
        identification = get(identification_rank),
        parent = ifelse(is.na(parent_rank), NA, get(parent_rank)),
        grandparent = ifelse(is.na(grandparent_rank), NA, get(grandparent_rank))
    ) %>%  
    dplyr::ungroup() 

# join to ncbi_lineageparents based on the three identification taxa and their ranks (replace BOLD with NCBI lineage)
bold_ncbi_joined <- 
    bold_id_tibble %>%
    dplyr::left_join(
        ., 
        ncbi_lineageparents, 
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
        nuc,
        identification_method,
        insdc_acs,
        sequence_run_site
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
        nuc,
        identification_method,
        insdc_acs,
        sequence_run_site
    ) %>%
    dplyr::mutate(
        taxid = stringr::str_replace(taxid, "^", "NCBI:"),
        across(kingdom:species, .fns = ~replace(., is.na(.), "Unclassified"))
    )

# combine matched and unmatched
bold_seqs <- 
    dplyr::bind_rows(bold_matched_seqs, bold_unmatched_seqs)

## filter and process the sequence records before exporting
bold_seqs_processed <- 
	bold_seqs %>% 
	# make sure all sequences are upper-case and convert I to N
	dplyr::mutate(
        nuc = toupper(nuc) %>% stringr::str_replace_all(., "I", "N")
    ) %>%
	# conditionally replace 'placeholder' species names (eg. "Genus sp. XYZ") with "Unclassified"
    {
        if (placeholder_as_unclassified) {
            dplyr::mutate(
                .,
                species = dplyr::case_when(
                    stringr::str_detect(species, "\\.") ~ "Unclassified", # if species name has ".", treated as unclassified
                    .default = species 
                )
            )
        } else { . }
    } %>%
    # conditionally replace digit-containing taxon names with "Unclassified"
    {
        if (digits_as_unclassified) {
            dplyr::mutate(
                .,
                dplyr::across(
                    kingdom:species,
                    ~dplyr::case_when(
                        stringr::str_detect(., "\\d") ~ "Unclassified", # if any taxonomic rank has at least one digit, treated as unclassified
                        .default = . 
                    )
                )
                
            )
        } else { . }
    } %>%
	# remove sequences mined from GenBank if GenBank is being used
    {	if (use_genbank){
    		dplyr::filter(., !stringr::str_detect(sequence_run_site, "GenBank|NCBI"))
    	} else { . }
    } %>%
	# conditionally filter out explicitly BOLD-classified sequences
    {   if( bold_idmethod_filter ) {
        	dplyr::filter(., !stringr::str_detect(identification_method,            "BIN|(B|b)arcode|DNA|BOLD|(T|t)ree|(S|s)equence"))
    	} else { . } 
    } %>%
	# categorise sequences into those that have explicit GenBank duplicates and those that don't
	dplyr::mutate(gb_dup = insdc_acs %in% gb_acc) %>%
	dplyr::select(-identification_method, -sequence_run_site)

## create FASTA format from sequence tibble
bold_seqs_prefasta <- 
    bold_seqs_processed %>%
    tidyr::unite("id", c(seqid, taxid), sep = "|") %>% # combine ids into a single column
    tidyr::unite("ranks", kingdom:species, sep = ";") %>% # combine ranks into a single column
    dplyr::mutate(
        ranks = stringr::str_replace_all(ranks, " +", " ") # replace runs of spaces with a single space (to allow HMMER output parsing later)
    ) %>%
    tidyr::unite("header", id:ranks, sep = ";") 
  
# DSS object of sequences (including GB-duplicates -- these get processed later)
seqs_out <- 
	bold_seqs_prefasta %>%
	dplyr::select(-gb_dup, -insdc_acs) %>%
	tibble::deframe() %>%
	Biostrings::DNAStringSet(.)

# table of GB-duplicate sequence names with their GB accessions
seqs_dup <- 
	bold_seqs_prefasta %>%
	dplyr::filter(gb_dup) %>%
	dplyr::select(-gb_dup, -nuc) %>%
	dplyr::rename(name = header, genbank_acc = insdc_acs)

### save outputs
# tibble of changes to BOLD sequence taxonomy by matching to NCBI synonyms
readr::write_csv(bold_ncbi_synchanges, "synchanges.csv")
# tibble of taxa (from sequences) that can be found in BOLD and NCBI 
readr::write_csv(matching_taxids, "matching_taxids.csv")
# write FASTA of non-GB-duplicate sequences 
Biostrings::writeXStringSet(seqs_out, "bold_seqs.fasta")
# write table of GB-duplicate sequences 
readr::write_csv(seqs_dup, "bold_gp_dups.csv")

# clean up environment before saving
rm(seqs_dup)
rm(seqs_out)
rm(bold_seqs_processed)
rm(bold_seqs_prefasta)
rm(bold_seqs)
rm(bold_matched_seqs)
rm(bold_ncbi_unmatched)
rm(bold_ncbi_joined)
rm(bold_id_tibble)
rm(bold_db_targets_syns)
rm(bold_db_targets)
gc()