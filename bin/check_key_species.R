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
    "ggplot2",
    "httr",
    "kmer",
    "magrittr",
    "methods",
    "openssl",
    "parallel",
    "patchwork",
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
    "sf_meta_file",
    "key_species_list"
    )
lapply(nf_vars, nf_var_check)

### process variables 
sf_meta <- readr::read_csv(sf_meta_file)

key_species <- 
    readr::read_lines(key_species_list) %>% 
    stringr::str_trim() %>% 
    as.list()

# define fate columns 
cols_fates <- c(
    filter_unclassified = 0,
    filter_phmm_full = 0,
    filter_phmm_trimmed = 0, 
    filter_duplicates = 0, 
    filter_ambiguous = 0, 
    filter_tax_outliers = 0, 
    filter_seq_outliers = 0, 
    filter_redundant = 0, 
    final_database = 0
)

### run code

# function that find sequences per species
find_key_sequences <- function(species_id, seq_tibble, add_root) {
    
    # check if species_id is a taxid or a name
    if (species_id %>% stringr::str_detect("^\\d+$")){
        # taxid
        species_id_valid <- species_id
        # get all sequences matching species
        species_seq_tibble <- 
            seq_tibble %>%
            dplyr::filter(stringr::str_detect(taxid, paste0("NCBI:",species_id_valid,"$")))
    } else {
        # string
        # convert species name into 'valid' internal format 
        species_id_valid <- 
        species_id %>% 
        stringr::str_replace_all(., " +", " ") # replace two or more spaces
        
        # get all sequences matching species
        species_seq_tibble <- 
            seq_tibble %>%
            dplyr::filter(species == species_id_valid)
    }
    # get vector of species sequences in the final sequences
    final_sequences <- 
        species_seq_tibble %>%
        dplyr::filter(fate == "final_database") %>%
        dplyr::pull(name)
    
    # add Root to lineage string if parameter is true
    if(add_root){ 
        final_sequences <- 
        final_sequences %>%
        stringr::str_replace(., ";", ";Root;")
    }
    
    if (rlang::is_empty(final_sequences)){
        final_sequences <- NULL
    }
    
    # make summary tibble
    if(is.null(final_sequences)){
        species_summary_tibble <- 
        tibble::tibble_row(
            kingdom = NA,
            phylum = NA,
            class = NA,
            order = NA, 
            family = NA,
            genus = NA, 
            species = species_id_valid,
            filter_unclassified = 0,
            filter_phmm_full = 0,
            filter_phmm_trimmed = 0, 
            filter_duplicates = 0, 
            filter_ambiguous = 0, 
            filter_tax_outliers = 0, 
            filter_seq_outliers = 0, 
            filter_redundant = 0, 
            final_database = 0
        )
        
    } else {
        species_summary_tibble <- 
        species_seq_tibble %>%
        dplyr::group_by(kingdom, phylum, class, order, family, genus, species, fate) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider(names_from = fate, values_from = n, values_fill = 0) %>%
        tibble::add_column(!!!cols_fates[!names(cols_fates) %in% names(.)]) %>%
        dplyr::relocate(
            filter_unclassified,
            filter_phmm_full,
            filter_phmm_trimmed,
            filter_duplicates,
            filter_ambiguous,
            filter_tax_outliers,
            filter_seq_outliers,
            filter_redundant,
            final_database,
            .after = species
        )
    }
    
    # return list 
    return(
        list(
        species_id, 
        species_seq_tibble, 
        species_summary_tibble,
        final_sequences
        )
    )
}

# get outputs per key species
key_species_output <- 
    lapply(
        key_species,
        find_key_sequences, 
        seq_tibble = sf_meta, 
        add_root = TRUE
    )

# function that extracts nth sub-element of a each element of a list
extract_subelement <- function(lst, n){
    sapply(lst, `[`, n)
}

# table of all sequences assigned to key species
key_sequences <- 
    extract_subelement(key_species_output, 2) %>% 
    dplyr::bind_rows()

readr::write_csv(key_sequences, "key_sequences.csv")

# summary of number of sequences per fate for each key species
key_species_summary <- 
  extract_subelement(key_species_output, 3) %>% 
  dplyr::bind_rows()

readr::write_csv(key_species_summary, "key_species_summary.csv")

# saving final database sequences per key species 
extract_subelement(key_species_output, 4) %>% 
    .[!sapply(.,is.null)] %>% # remove null elements
    lapply(
        .,
        function(x){
            # get name of species
            spp_name <- x[1] %>% stringr::str_extract(., "(?<=;)[^;]+?$")

            # save vector as text file
            readr::write_lines(x, paste0(spp_name,".seqnames.txt"))

            print(paste0("Wrote sequence names for ",spp_name," to file"))

            return(NULL)
        }
    )
  
  

