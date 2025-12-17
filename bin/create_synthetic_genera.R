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
    "fasta_files",
    "clusters_tsv"
)
lapply(nf_vars, nf_var_check)

### process variables 

clusters <- readr::read_tsv(clusters_tsv, col_names = c("rep","name"), show_col_types = FALSE)

# read in list of sequences
fasta_list <-
    stringr::str_extract_all(fasta_files, pattern = "[^\\s,\\[\\]]+") %>% unlist()

seqs_list <- 
    lapply(
        fasta_list, 
        Biostrings::readDNAStringSet
    ) %>%
    lapply(
        .,
        function(x){
            # within each cluster, extract the genus lineage and rename
            lin <- 
                names(x) %>%
                tibble::as_tibble_col(column_name = "name") %>%
                dplyr::mutate(genus_old = stringr::str_extract(name, "(?<=;).+$") %>% stringr::str_extract(., ".+(?=;)")) %>%
                # add cluster reps
                dplyr::left_join(., clusters, by = "name") %>%
                # add unique id for each cluster to the genus name
                dplyr::group_by(rep) %>%
                dplyr::mutate(
                    genus_new = paste0(genus_old, dplyr::cur_group_id())
                ) %>%
                dplyr::ungroup() %>%
                # create a new name with new genus
                dplyr::mutate(
                    name_new = stringr::str_replace(name, genus_old, genus_new)
                )
            # rename sequences
            names(x) <- lin$name_new
            
            return(x)
        }
    )

# collect all sequences into a single file
seqs <- seqs_list %>% do.call(c, .)

# write sequences to file
Biostrings::writeXStringSet(seqs, "seqs_out.fasta")