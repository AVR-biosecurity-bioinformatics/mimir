#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dplyr",
    "magrittr",
    "stringr",
    "tidyr",
    # "taxreturn",
    # "XVector",
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
    "genetic_code"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read in seqs from file
seqs <- readRDS(seqs_file)

### run code

## filter out sequences with stop codons
seqs_filtered <- 
    codon_filter(
        x = seqs, 
        genetic_code = genetic_code, 
        tryrc = TRUE, 
        resolve_draws = "majority"
    )

# save filtered sequences as .rds file
saveRDS(seqs_filtered, paste0(taxon,"_",type,"_filter_stop.rds"))

