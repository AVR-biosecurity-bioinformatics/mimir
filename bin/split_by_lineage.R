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
    "split_rank_input"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# convert fasta to DNAbin
seqs <- ape::read.FASTA(fasta_file)

# parse split_rank parameter
split_rank <- split_rank_input

# get index of split_rank for grouping below
rank_cols <- c("kingdom","phylum", "class", "order", "family", "genus", "species")
split_rank_index <- match(split_rank, rank_cols) 

# get the taxonomic lineage of each sequence as a tibble/df
lineage <- 
    names(seqs) %>%
    tibble::as_tibble(column_name = "value") %>%
    tidyr::separate(col=value, into=c("id", "lineage_string"), sep=";", extra="merge") %>%
    tidyr::separate(col=lineage_string, into=c("kingdom","phylum", "class", "order", "family", "genus", "species"), sep=";", extra="merge") %>%
    tidyr::separate(col=id, into=c("seqid", "taxid"), sep="\\|", extra="merge") 

# want to get all groups for each taxon combination down to the `split_rank` level
# eg. if split_rank is genus, want all lineages at the genus level 
lineage_grouped <- 
  lineage %>% 
  dplyr::group_by(dplyr::across(3:(split_rank_index+2))) # use index of split_rank to choose grouping columns

rm(lineage) # remove unneeded object

# split into list of tibbles by lineage ending at split_rank
lineage_split <- 
  lineage_grouped %>%
  dplyr::group_split()
  
# regenerate header structure
split_list <- 
  lapply(lineage_split, tidyr::unite, col = "ID", seqid:taxid, sep = "|") %>%
  lapply(., tidyr::unite, col = "lineage", kingdom:species, sep = ";") %>%
  lapply(., tidyr::unite, col = "header", ID:lineage, sep = ";") %>%
  lapply(., dplyr::pull, "header")

rm(lineage_split) # remove unneeded object

# name list  using group keys
names(split_list) <- 
  lineage_grouped %>% 
  dplyr::group_keys() %>%
  tidyr::unite(col = "name", sep = "_") %>%
  dplyr::pull(name)

rm(lineage_grouped) # remove unneeded object

# function to subset DNAbin object based on named list of sequence names, saving as fasta
fasta_from_list <- function(named_list, DNAbin) {
  # vector of sequence names in input list
  seq_names <- named_list %>% purrr::as_vector() %>% unname()
  # subset DNAbin
  retained <- DNAbin[names(DNAbin) %in% seq_names]
  return(retained)
}

# create list of sequence subsets
seqs_subsets <-
  sapply(
    split_list, 
    fasta_from_list, 
    DNAbin = seqs, 
    simplify = FALSE, 
    USE.NAMES = TRUE
    )

rm(split_list) # remove unneeded object

# save fasta of each subset
sapply(
    names(seqs_subsets),
    function (name) {
        ape::write.FASTA(
            x = seqs_subsets[[name]],
            file = paste0(name,".grouped.fasta")
        )
    },
    simplify = FALSE, 
    USE.NAMES = TRUE
)