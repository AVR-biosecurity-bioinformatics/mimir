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
    "taxon",
    "taxon_rank",
    "entrez_key",
    "ncbi_synonyms_file",
    "ncbi_gencodes_file"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# set API key if provided
if ( entrez_key != "no_key" ){
    Sys.setenv(ENTREZ_KEY = entrez_key)
}

# process taxize 'rank_ref' df into vector of valid ranks
valid_ranks <- 
    rank_ref$ranks %>%
    stringr::str_split(",") %>%
    unlist()

# taxon_rank
if ( taxon_rank == "no_ranks" ){
    rank_query <- NULL
} else {
    # check rank is valid
    if ( taxon_rank %in% valid_ranks ) {
        rank_query <- taxon_rank
    } else {
        stop ( paste0( "ERROR: '--target_ranks value' ('",taxon_rank,"') is not a valid taxonomic rank.\n\tSelect one of: ",stringr::str_c(valid_ranks, collapse = ', ')))
    }
}

# import ncbi synonyms tibble
ncbi_synonyms <- readRDS(ncbi_synonyms_file)

# import ncbi synonyms tibble
ncbi_gencodes <- readRDS(ncbi_gencodes_file)

# sorted NCBI ranks from Supplemental Table 3 in Schoch et al 2020 (https://doi.org/10.1093%2Fdatabase%2Fbaaa062)
sorted_ranks <- c(
    "superkingdom",
    "kingdom",
    "subkingdom",
    "superphylum",
    "phylum",
    "subphylum",
    "infraphylum",
    "superclass",
    "class",
    "subclass",
    "infraclass",
    "cohort",
    "subcohort",
    "superorder",
    "order",
    "suborder",
    "infraorder",
    "parvorder",
    "superfamily",
    "family",
    "subfamily",
    "tribe",
    "subtribe",
    "genus",
    "subgenus",
    "section",
    "subsection",
    "series",
    "subseries",
    "species group",
    "species subgroup",
    "species",
    "forma specialis",
    "subspecies",
    "varietas",
    "subvariety",
    "forma",
    "serogroup",
    "serotype",
    "strain",
    "isolate"
)

# BOLD valid taxonomic ranks
bold_validranks <- c(
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "subfamily",
    "tribe",
    "genus",
    "species"
)

### run code

## check if taxon is NCBI UID (string of digits) 
if ( stringr::str_detect(taxon, "^\\d+$") ) { # if string of only digits...
    # convert id to numeric
    ncbi_id <- as.numeric(taxon)
    # get character taxon name from id
    name_from_ncbi_id <- taxize::id2name(ncbi_id, db = "ncbi")
    # check id is valid
    if ( identical(name_from_ncbi_id[[1]]$name, character(0) ) ){
        stop (paste0( "ERROR: --target_taxa value '",ncbi_id,"' is not a valid NCBI taxid (or there is an issue with NCBI)." ))
    } else {
        taxon_name <- name_from_ncbi_id[[1]]$name
    }
    # check rank of id is that same as that supplied to pipeline
    if ( taxon_rank != "no_ranks" ){
        if ( !identical(name_from_ncbi_id[[1]]$rank, taxon_rank)) {
            stop( paste0("ERROR: Taxonomic rank of supplied taxon ID ('",taxon,"') is different to '--target_ranks' value ('",taxon_rank,"').") ) 
        }
    } else {
        # get rank from NCBI query
        taxon_rank <- name_from_ncbi_id[[1]]$rank
    }
} else { # otherwise assume taxon supplied is a name
    # convert taxon name to ncbi uid
    ncbi_id <- 
        taxize::get_uid(
            sci_com = taxon,
            modifier = "Scientific Name",
            rank_query = rank_query, # narrow down to this rank, if supplied
            ask = FALSE # make function return "NA" if there is ambiguity
        )[1] 
    # check supplied taxon name returns valid uid
    if (is.na(ncbi_id)){
        stop ( paste0("ERROR: '--target_taxa' value ('",taxon,"') cannot be found in the NCBI Taxonomy database (or there is an issue with NCBI)."))
    } else {
        ncbi_id <- as.numeric(ncbi_id) # force numeric class
    }
    # get rank if not supplied
    if ( taxon_rank == "no_ranks" ){
        taxon_rank <- 
            taxize::tax_rank(
                sci_id = taxon,
                db = "ncbi"
            )[[1]]
    }
    # save taxon name variable
    taxon_name <- taxon
}

## chunk to downstream valid ranks when target taxon is a rank that is outside of the BOLD ranks 
if (!taxon_rank %in% bold_validranks){
    ## find highest "valid" BOLD rank downstream of target taxon 
    # get index of target rank in sorted_ranks
    taxon_rank_index <- match(taxon_rank, sorted_ranks) 
    # get indices of bold_validranks in sorted_ranks
    bold_validranks_index <- match(bold_validranks, sorted_ranks)
    # get name of BOLD rank
    downstream_valid_rank <- bold_validranks[min(which(bold_validranks_index > taxon_rank_index))]

    ## split into lower ranks in NCBI, then search for those in BOLD as chunks
    # get vector of downstream NCBI ids to search BOLD 
    rank_chunks <- 
        taxize::ncbi_downstream(
            id =  ncbi_id,
            downto = downstream_valid_rank
        )
    if (length(rank_chunks) == 0){stop("Retrieving downstream NCBI ids failed; check ENTREZ key (or there is an issue with NCBI).")}
    # create search pattern from NCBI ids
    rank_chunk_pattern <- 
        rank_chunks$childtaxa_id %>%
        stringr::str_replace(pattern = "^(.*?)$", replacement = "^\\1$")
    # search NCBI synonyms list to get vector of possible alternative matches at the correct rank
    rank_chunk_synonyms <- 
        ncbi_synonyms %>%
        dplyr::filter(
            stringr::str_detect(tax_id, stringr::str_c(rank_chunk_pattern, collapse = "|")) & 
            rank == downstream_valid_rank
        ) %>%
        dplyr::pull(name_txt)
    # combine NCBI names and synonyms into one vector
    rank_chunks_complete <- c(rank_chunks$childtaxa_name, rank_chunk_synonyms)
    # search BOLD with vector of downstream taxa from NCBI
    rank_chunks_boldids <- 
        taxize::get_boldid_(
            sci = rank_chunks_complete, 
            rank = downstream_valid_rank
        )
    # combine data frames, removing rows of taxa that aren't in BOLD (taxid is NA)
    rank_chunks_df <- 
        dplyr::bind_rows(rank_chunks_boldids) %>% 
        dplyr::filter(!is.na(taxid))
    # remove rows that have their parent taxon in the list of taxa
    rank_chunks_df_final <- 
        rank_chunks_df[is.na(with(rank_chunks_df, match(parentname, taxon))),] 
    # get BOLD taxon names
    bold_ids <- rank_chunks_df_final$taxid
    # get BOLD taxon ids
    bold_names <- rank_chunks_df_final$taxon
    # get BOLD taxon rank
    bold_rank <- downstream_valid_rank

} else { ## if target rank is in BOLD ranks...
    ## convert taxon name to BOLD ID
    bold_ids <- 
        taxize::get_boldid(
            sci = taxon_name,
            messages = FALSE,
            fuzzy = FALSE,
            rank = rank_query
        )[1] %>%
            as.integer()

    ## convert BOLD ID to BOLD taxon name
    bold_names <-
        taxize::id2name(
            id = bold_ids, 
            db = "bold"
        )[[1]]$name

    # check BOLD and NCBI taxon names are the same
    if (!identical(taxon_name, bold_names)) {
        stop ( paste0("ERROR: Taxonomic names in NCBI ('",taxon_name,"') and BOLD ('",bold_names,"') are not identical!"))
    }

    # get BOLD taxon rank 
    bold_rank <- taxon_rank
}

### get genetic codes for taxon (TODO: add support for lists of taxa, using taxize::lowest_common)
taxon_gencodes <- 
    ncbi_gencodes %>%
    dplyr::filter(tax_id == ncbi_id) %>%
    dplyr::select(tax_id, tax_name, gencode, mtgencode, plgencode, hdgencode)

### save output files
# NCBI name
write(taxon_name, paste0(taxon,"_name.txt"), ncolumns = 1, sep = "\n")
# NCBI ID
write(ncbi_id, paste0(taxon,"_ncbi_id.txt"), ncolumns = 1, sep = "\n")
# BOLD name(s)
write(bold_names, paste0(taxon, "_bold_names.txt"), ncolumns = 1, sep = "\n")
# BOLD ID(s)
write(bold_ids, paste0(taxon,"_bold_ids.txt"), ncolumns = 1, sep = "\n")
# BOLD taxon rank
write(bold_rank, paste0(taxon, "_bold_rank.txt"), ncolumns = 1, sep = "\n")
# gencodes for rank
readr::write_delim(taxon_gencodes, paste0(taxon, "_gencodes.csv"), delim = ",")