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
    "phmm_model_file",
    "coding",
    "task_index"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read in phmm model from file
phmm_model <- readRDS(phmm_model_file)

# convert fasta to DNAbin
seqs <- ape::read.FASTA(fasta_file, type = "DNA")

## parameter parsing
min_score <- as.numeric(params.phmm_min_score)

min_length <- as.numeric(params.phmm_min_length)

if ( params.shave_to_phmm == "true" ) {
    shave <- TRUE
} else {
    shave <- FALSE
}

if ( coding == "true" ) { 
    check_frame <- TRUE 
} else {
    check_frame <- FALSE
}

if ( params.remove_ambiguous == "true" ) {
    max_N <- 0
} else {
    max_N <- Inf 
}

### run code

## filter using PHMM model
seqs_filtered <- 
    map_to_model(
        x = seqs, # sequences, as a DNAbin or DNAStringset object
        model = phmm_model, # PHMM, as a PHMM object
        min_score = min_score, # minimum specificity of match for retention (see ?aphid::Viterbi)
        min_length = min_length, # minimum length of match for retention 
        max_N = max_N, # max ambiguous bases to allow 
        max_gap = Inf, # max gaps to allow
        shave = shave, # whether to remove bases outside match to PHMM
        check_frame = check_frame, # check if indels are in multiples of three
        kmer_threshold = 0.5, # k-mer distance allowed for retention for alignment step
        k = 5, # k-mer size for clustering
        multithread = FALSE, 
        quiet = FALSE
    )

# write fasta
if ( !is.null(seqs_filtered) && length(seqs_filtered) > 0 ){
    write_fasta(
        seqs_filtered, 
        file = paste0("filter_phmm.",task_index,".fasta"), 
        compress = FALSE
        )
} else {
    file.create(paste0("filter_phmm.",task_index,".fasta"))
}