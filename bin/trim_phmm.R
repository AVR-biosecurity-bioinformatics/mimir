#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dplyr",
    "magrittr",
    "stringr",
    "tidyr",
    # "taxreturn",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "params_dict",
    "phmm_model_file",
    "primer_fwd",
    "primer_rev",
    "remove_primers"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read in phmm model from file
phmm_model <- readRDS(phmm_model_file)

# parse remove_primers
if ( remove_primers == "true" ) {
    trimprimers <- TRUE
} else if ( remove_primers == "false" || remove_primers == "null" ) {
    trimprimers <- FALSE
} else {
    stop(paste0("ERROR: params.remove_primers is '",remove_primers,"' but must be 'true', 'false' or 'null'."))
}

### run code

## trim PHMM model
phmm_model_trimmed <- 
    subset_model(
        x = phmm_model, 
        primers = c(primer_fwd, primer_rev), 
        trimprimers = trimprimers
    )

# save trimmed model as .rds file
saveRDS(phmm_model_trimmed, "phmm_model_trimmed.rds")

