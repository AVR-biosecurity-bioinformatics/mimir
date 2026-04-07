#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "Biostrings",
    "DECIPHER",
    "IRanges",
    "R.utils",
    "RCurl",
    "ape",
    # "aphid",
    "bold",
    "data.table",
    "data.tree",
    "dplyr",
    # "entropy",
    # "fs",
    # "furrr",
    # "future",
    "httr",
    "kmer",
    "magrittr",
    "methods",
    "openssl",
    "parallel",
    # "phytools",
    "purrr",
    "readr",
    # "rentrez",
    # "rvest",
    "stats",
    "stringr",
    # "taxize",
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
    "fasta_file"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read in sequences
seqs <- ape::read.FASTA(fasta_file)

# length of each sequence
seq_lengths <- 
    seqs %>% 
    DNAbin2DNAstringset() %>% 
    width()

# get lineage info with sequence lengths
seqs_info <-  
    seqs %>%
    names() %>%
    tibble::as_tibble(column_name = "value") %>%
    tidyr::separate(col=value, into=c("id", "lineage_string"), sep=";", extra="merge") %>%
    tidyr::separate(col=lineage_string, into = c("kingdom","phylum", "class", "order", "family", "genus", "species"), sep=";", extra="merge") %>%
    tidyr::separate(col=id, into=c("seqid", "taxid"), sep="\\|", extra="merge") %>%
    dplyr::mutate(length = seq_lengths)

# get number of families in dataset
n_families <- 
    seqs_info %>% 
    dplyr::pull(family) %>% 
    unique() %>% 
    length()

# calculate maximum number of sequences to extract from each family to get core to ~1000 sequences
n_extract <- ceiling(1000/n_families)

# get median length of sequences
length_median <- 
    seqs_info %>% 
    dplyr::pull(length) %>% 
    median()

# select longest sequence from each family
seqs_core_ids <- 
    seqs_info %>%
    dplyr::group_by(family) %>%
    dplyr::arrange(desc(length)) %>%
    dplyr::slice(1:n_extract) %>%
    dplyr::ungroup() %>%
    # remove sequences below the median length
    dplyr::filter(length >= length_median) %>%
    # unite seqid and taxid, then pull
    tidyr::unite(col = "id", seqid:taxid, sep = "|") %>%
    dplyr::pull(id)

# select 200 longest sequences overall
seqs_longest_ids <- 
    seqs_info %>%
    dplyr::arrange(desc(length)) %>%
    dplyr::slice(1:200) %>%
    tidyr::unite(col = "id", seqid:taxid, sep = "|") %>%
    dplyr::pull(id)

# make search pattern from ids
seqs_core_pattern <- 
    c(seqs_core_ids, seqs_longest_ids) %>% # combine longest and longest per family
    stringr::str_unique() %>% 
    stringr::str_replace_all("\\|", "\\\\|") %>% # escape |
    stringr::str_flatten(collapse = "|")

# extract core sequences
seqs_core <- seqs[names(seqs) %>% stringr::str_starts(., seqs_core_pattern)]

# extract other sequences
seqs_other <- seqs[!names(seqs) %>% stringr::str_starts(., seqs_core_pattern)]

# write sequences to file
ape::write.FASTA(seqs_core, "core.fasta")

ape::write.FASTA(seqs_other, "other.fasta")