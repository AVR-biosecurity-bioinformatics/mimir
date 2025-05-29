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
    "db_tsv_file",
    "db_meta_file",
    "bold_tibble_file",
    "marker",
    "bold_idmethod_filter"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read in bold db
bold_db <- 
    readr::read_delim(
        db_tsv_file, 
        delim = "\t"
    )

# read in bold names and rank as tibble (and deduplicate)
bold_tibble <- readr::read_csv(bold_tibble_file) %>% dplyr::distinct()

# min and max sequence lengths
min_length <- as.integer(params.min_length)
max_length <- as.integer(params.max_length)

if ( bold_idmethod_filter == "true" ){
    bold_idmethod_filter <- TRUE
} else {
    bold_idmethod_filter <- FALSE
}

### run code

# get subset tibbles that match taxon in each row of the bold_tibble
list_out <- purrr::imap(
    .x = bold_tibble %>% deframe() %>% as.list(),
    .f = \(x, idx, db = bold_db){
        output <- 
            db %>%
            dplyr::filter(get({{x}}) == {{idx}})
        return(output)
    }
)

# subset BOLD db
bold_db_targets <- 
    list_out %>%
    dplyr::bind_rows() %>%
    # get only marker of interest
    dplyr::filter(marker_code == marker) %>%
    # remove rows where there is no sequence
    dplyr::filter(!is.na(nuc)) %>%
    # remove sequences mined from GenBank
    dplyr::filter(!stringr::str_detect(sequence_run_site, "GenBank|NCBI")) %>% 
    # conditionally filter out explicitly BOLD-classified sequences
    {   if( bold_idmethod_filter )
            dplyr::filter(., !stringr::str_detect( # note the '.' is essential here
                identification_method, 
                "BIN|(B|b)arcode|DNA|BOLD|(T|t)ree|(S|s)equence"
                )
            )
        else 
            . 
    } %>%
    # dealign nucleotides
    dplyr::mutate(nuc = stringr::str_remove_all(nuc, "\\-")) %>%
    # filter nucleotide length
    dplyr::mutate(nuc_basecount = as.integer(nuc_basecount)) %>%
    dplyr::filter(nuc_basecount >= min_length, nuc_basecount <= max_length) %>%
    # select only columns of interest
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
        nuc
    ) %>%
    # one last deduplication just in case
    dplyr::distinct()

## save subset database as .rds (if data still remains after filtering)
if ( nrow(bold_db_targets) > 0 ){
    saveRDS(bold_db_targets, paste0("bold_db_targets.rds")) ### TODO: make this naming work with multiple taxa as pipeline input
} else {
    message("No sequences matching criteria found -- not saving output .rds file.")
}

