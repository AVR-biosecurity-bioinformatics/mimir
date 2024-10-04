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
    "entrez_key",
    "db_path"
    )
lapply(nf_vars, nf_var_check)

### process variables 

chunk_rank <- params.chunk_rank

# set API key if provided
if ( entrez_key != "no_key" ){
    Sys.setenv(ENTREZ_KEY = entrez_key)
}

## convert taxon name to uid, if necessary
if ( stringr::str_detect(taxon, "^\\d+$") ) { # if already a number, leave as is
    # ensure numeric class
    id <- as.numeric(taxon)
    # get character taxon name from id
    taxon_name <- taxize::id2name(id, db = "ncbi")[[1]] %>% pull(name)
} else {
    # convert taxon name to uid
    id <- 
        taxize::get_uid(
            sci_com = taxon,
            modifier = "Scientific Name"
        )[1]
    # save taxon name variable
    taxon_name <- taxon
}

if ( is.na(id) ){
    stop( paste0("ERROR: Target taxon '",taxon,"' could not be unambiguously converted to UID/taxid.\n\tConsider supplying UID directly."))
}

# import 'taxidlineage.dmp' from taxdump
taxidlineage.dmp <- 
    readr::read_tsv(
        paste0(db_path,"/taxidlineage.dmp"),
        col_names = c("tax_id", "lineage"),
        col_types = ("i-c-")
    )

# import 'nodes.dmp' from taxdump
nodes.dmp <- 
    readr::read_tsv(
        paste0(db_path,"/nodes.dmp"),
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

# check chunk_rank is a valid rank
all_ranks <- 
    nodes.dmp %>%
    dplyr::pull(rank) %>% 
    unique()

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


if ( !chunk_rank %in% all_ranks ){
    stop ( paste0("ERROR: params.chunk_rank ('",chunk_rank"') is not a valid NCBI taxonomic rank!\n\tSelect one of: ",stringr::str_c(all_ranks, collapse = ', ')))
}

### run code


## check if target taxon rank is higher than the chunk rank
# get classification of target rank
taxon_classification <- 
    taxize::classification(
        sci_id = id,
        db = "ncbi"
    )
# get actual rank and upstream ranks of target taxon
taxon_ranks <- taxon_classification[[1]] %>% pull(rank)

# subset taxidlineage.dmp to taxid that contain target taxon in their upstream lineage
downstream_taxids <- 
    taxidlineage.dmp %>%
    dplyr::filter(
        stringr::str_detect(
            lineage,
            paste0("( ",id," )|( ",id,"$)")
        )
    )

# get all nodes that are of chunk_rank
chunk_rank_nodes <- 
    nodes.dmp %>%
        dplyr::filter(rank == chunk_rank) %>%
        dplyr::pull(tax_id)

# get nodes within target that are of chunk_rank
target_chunks <- 
    downstream_taxids %>%
        dplyr::filter(tax_id %in% chunk_rank_nodes) %>%
        dplyr::pull(tax_id)

## find downstream nodes that don't have assigned taxonomy at chunk_rank (to fetch otherwise missing sequences)
# create search patterns from target chunks
target_chunk_pattern <- # match internal taxids in lineage string
    stringr::str_c(target_chunks, collapse = " | ") %>%
        stringr::str_replace_all("^|$"," ") 

target_chunk_pattern_endmatch <- # match end taxids in lineage string
    target_chunk_pattern %>%
    stringr::str_replace_all(" \\|","$|") %>%
    stringr::str_replace(" $","$")

# get nodes that don't have target_chunks in lineage
nontarget_nodes <-
    downstream_taxids %>%
        dplyr::filter(!tax_id %in% target_chunks) %>% # remove taxids that match wanted chunks
        dplyr::filter( # remove taxids that contain target_chunks in their upstream lineage
            !stringr::str_detect(
                lineage, 
                target_chunk_pattern
            )
        ) %>%
        dplyr::filter( # remove taxids that contain target_chunks in their upstream lineage (endmatching)
            !stringr::str_detect(
                lineage, 
                target_chunk_pattern_endmatch
            )
        )

# get index of chunk_rank in sorted_ranks vector
chunk_rank_index <- match(chunk_rank, sorted_ranks)

# join nodes.dmp data to nontarget_nodes, get rank index, and remove nodes with rank higher than chunk_rank
nontarget_nodes %>%
    dplyr::left_join(., nodes.dmp, by = "tax_id") %>%
    dplyr::mutate(
        sorted_index = match(rank, sorted_ranks) # add index of rank to tibble
    ) %>%
    dplyr::filter(sorted_index > 20)


# ## if chunk rank is lower than taxon rank, chunk taxon
# if ( !is.element(params.chunk_rank, taxon_ranks ) ){
#     ## Fetch downstream taxa at a given rank
#     taxon_chunked <- 
#         taxize::ncbi_downstream(
#             id = id, 
#             downto = params.chunk_rank
#         )

#     # subset to just a vector of names
#     chunk_vector <- taxon_chunked$childtaxa_id    
# } else { 
#     ## if chunk rank same or higher than taxon rank, just use taxon as the output
#     chunk_vector <- taxon_name
# } 

# check that a vector was produced
if ( is.null(chunk_vector) ){
    stop (paste0("ERROR: No taxon list produced from taxid = '",id,"' and chunk_rank = '",params.chunk_rank,"'.\n\tCheck target taxon is not at or below the chunk rank."))
}

# save vector as text file
write(chunk_vector, "tax_list.txt")
