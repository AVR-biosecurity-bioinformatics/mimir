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
    "fasta_file",
    "subsample_seeds",
    "subsample_size"
    )
lapply(nf_vars, nf_var_check)

### process variables 

seeds <- readr::read_lines(subsample_seeds) %>% as.integer()
subsample_size <- as.numeric(subsample_size)

# read fasta from file
seqs <- Biostrings::readDNAStringSet(fasta_file)


### run code

# handle when there are fewer or the same number of sequences as the subsample size
if (length(seqs) <= subsample_size){
  subsample_size <- length(seqs)
}

for (i in seeds){
    set.seed(i)
    sel_seqs <- character()
    # get species-level lineages for all sequences
    seq_lin <- 
        seqs %>%
        names() %>%
        tibble::as_tibble_col(column_name = "name") %>%
        dplyr::mutate(
            lineage = stringr::str_extract(name, "(?<=;).+$")
        ) 
    # get unique lineages
    lineages <- 
        seq_lin %>%
        dplyr::distinct(lineage) %>%
        dplyr::pull(lineage)
    while (length(sel_seqs) < subsample_size){
        # select a random lineage
        x <- base::sample(lineages, size = 1, replace = T)
        # set sample size depending on how many sequences remain
        if (subsample_size - length(sel_seqs) == 1){
            samp_size <- 1
        } else {
            samp_size <- 2
        }
        # get up to 2 sequences from that lineage
        sel <- 
            seq_lin %>%
            dplyr::filter(lineage == x) %>%
            dplyr::slice_sample(n = samp_size) %>%
            dplyr::pull(name)
        # remove sequences from tibble
        seq_lin <- 
            seq_lin %>%
            dplyr::filter(!name %in% sel)
        # recalculate lineages
        lineages <- 
            seq_lin %>%
            dplyr::distinct(lineage) %>%
            dplyr::pull(lineage)
        # add sequences to list
        sel_seqs <- c(sel_seqs, sel)
            
    }
    # get final sequence subsample
    seqs_sub <- seqs[names(seqs) %in% sel_seqs]
    # write sequences to file
    Biostrings::writeXStringSet(seqs_sub, stringr::str_glue("subsample{i}.fasta"))
}
