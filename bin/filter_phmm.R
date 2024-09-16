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
    "phmm_model_file"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read in phmm model from file
phmm_model <- readRDS(phmm_model_file)

# convert fasta to DNAbin
seqs <- 
    ape::read.FASTA(fasta_file, type = "DNA")

# get basename
file_basename_noext <- 
    basename(fasta_file) %>%
    stringr::str_remove("\\.fasta$")

### run code

## filter using PHMM model
seqs_filtered <- 
    map_to_model(
        x = seqs, # sequences, as a DNAbin or DNAStringset object
        model = phmm_model, # PHMM, as a PHMM object
        min_score = 100, # minimum specificity of match for retention (see ?aphid::Viterbi)
        min_length = 100, # minimum length of match for retention 
        max_N = Inf, # max ambiguous bases to allow 
        max_gap = Inf, # max gaps to allow
        shave = TRUE, # whether to remove bases outside match to PHMM
        check_frame = TRUE, # check if indels are in multiples of three
        kmer_threshold = 0.5, # k-mer distance allowed for retention for alignment step
        k = 5, # k-mer size for clustering
        multithread = FALSE, 
        quiet = FALSE
    )

# save filtered sequences as .rds file
if ( !is.null(seqs_filtered)) {
    saveRDS(seqs_filtered, paste0(file_basename_noext,"_filter_phmm.rds"))
}

# write fasta for debugging
if ( params.all_fasta == "true" && !is.null(seqs_filtered) ){
    write_fasta(
        seqs_filtered, 
        file = paste0(file_basename_noext,"_filter_phmm.fasta"), 
        compress = FALSE
        )
}