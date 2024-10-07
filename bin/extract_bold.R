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
    "bold_names_file",
    "bold_rank_file",
    "marker"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# get name of db_tsv chunk for output naming
chunk_index <- tools::file_path_sans_ext(basename(db_tsv_file)) 

# read in bold db
bold_db <- 
    readr::read_delim(
        db_tsv_file, 
        delim = "\t"
    )

# read in bold names and rank as vectors
bold_names <- readr::read_lines(bold_names_file)
bold_rank <- readr::read_lines(bold_rank_file)

### run code

# subset BOLD db
bold_db_targets <- 
    bold_db %>%
    # get only taxa of interest
    dplyr::filter(get({{bold_rank}}) == {{bold_names}}) %>%
    # get only marker of interest
    dplyr::filter(marker_code == marker) %>%
    # remove rows where there is no sequence
    dplyr::filter(!is.na(nuc)) %>%
    # remove sequences mined from GenBank
    dplyr::filter(!stringr::str_detect(sequence_run_site, "GenBank|NCBI")) %>% 
    # conditionally filter out explicitly BOLD-classified sequences
    {   if( params.bold_idmethod_filter == "true")
            dplyr::filter(., !stringr::str_detect( # note the '.' is essential here
                identification_method, 
                "BIN|(B|b)arcode|DNA|BOLD|(T|t)ree|(S|s)equence"
                )
            )
        else 
            . 
    } %>%
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
    )

## save subset database as .rds (if data still remains after filtering)
if ( nrow(bold_db_targets) > 0 ){
    saveRDS(bold_db_targets, paste0("bold_db_targets.",chunk_index,".rds")) ### TODO: make this naming work with multiple taxa as pipeline input
} else {
    message("No sequences matching criteria found -- not saving output .rds file.")
}

