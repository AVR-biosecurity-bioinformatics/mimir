#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dplyr",
    "magrittr",
    "stringr",
    "tidyr",
    "taxreturn",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "params_dict",
    "taxon",
    "type",
    "seqs_file",
    "phmm_model_file"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read in phmm model from file
phmm_model <- readRDS(phmm_model_file)

# read in seqs from file
seqs <- readRDS(seqs_file)

### run code

## filter using PHMM model
seqs_filtered <- 
    taxreturn::map_to_model(
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

# save trimmed model as .rds file
saveRDS(seqs_filtered, paste0(taxon,"_",type,"_filter_phmm.rds"))

