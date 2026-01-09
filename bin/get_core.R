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
    "clusters_tsv"
)
lapply(nf_vars, nf_var_check)

### process variables 

clusters <- readr::read_tsv(clusters_tsv, col_names = c("rep", "name"), show_col_types = F)
seqs <- Biostrings::readDNAStringSet(fasta_file)

### run code 

rep_names <- clusters$rep %>% unique()
n_seqs <- length(seqs)
# proportion of total sequences that are cluster reps
rep_prop <- length(rep_names) / n_seqs

# if reps are <5% of the sequences, top up to 5% of the total sequences with a sample of other sequences
if (rep_prop < 0.05){
	# remainder to get
	n_remainder <- ceiling((0.05 - rep_prop) * n_seqs)
	other_names <- names(seqs)[!names(seqs) %in% rep_names]
    # sample from other names
    remainder_names <- sample(other_names, n_remainder, replace = FALSE)
    # combine with cluster reps
    core_names <- c(rep_names, remainder_names)
} else {
	core_names <- rep_names
}

# get core and other
seqs_core <- seqs[names(seqs) %in% core_names]
seqs_other <- seqs[!names(seqs) %in% core_names]

# write to file
if (length(seqs_core) > 0){
	Biostrings::writeXStringSet(seqs_core, "seqs_core.fasta", width = 9999)
} else {
	stop("No core sequences identified -- please check input files")
}

if (length(seqs_other) > 0){
	Biostrings::writeXStringSet(seqs_other, "seqs_other.fasta", width = 9999)
} else {
	file.create("seqs_other.fasta")
}