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
    "ig_tibble_file",
    "cluster_reps_file"
)
lapply(nf_vars, nf_var_check)

### process variables 

seqs <- Biostrings::readDNAStringSet(fasta_file)

ig_tibble <- readr::read_csv(ig_tibble_file)

cluster_reps <- readr::read_lines(cluster_reps_file)

# get single rep from each consistent genus 
consistent_reps <- 
	ig_tibble %>%
	dplyr::filter(
		threshold == "genus_min", # keep genus_min rows
		name %in% names(seqs), # keep only sequences that passed genus and species filtering
		consistent == TRUE #keep only consistent genera
	) %>%
	dplyr::arrange(desc(n)) %>% # most frequent sequences on top
	dplyr::group_by(taxon) %>%
	dplyr::slice(1) %>%
	dplyr::ungroup() %>%
	dplyr::pull(name)

# get all sequences from inconsistent genera
inconsistent_names <- 
	ig_tibble %>%
	dplyr::filter(
		threshold == "genus_min", # keep genus_min rows
		name %in% names(seqs), # keep only sequences that passed genus and species filtering
		consistent == FALSE #keep only consistent genera
	) %>%
	dplyr::pull(name)

# get selected sequences
seqs_out <- seqs[names(seqs) %in% c(consistent_reps, inconsistent_names, cluster_reps)]

Biostrings::writeXStringSet(seqs_out, "reps.fasta")
