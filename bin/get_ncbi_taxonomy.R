#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dplyr",
    "magrittr",
    "purrr",
    "readr",
    "stringr",
    "tidyr",
    # "taxreturn",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "params_dict"
    )
lapply(nf_vars, nf_var_check)

### process variables 



### run code

# download ncbi taxonomy files
# db <- taxreturn::get_ncbi_taxonomy()
db <- get_ncbi_taxonomy()


# save db object
saveRDS(db, "ncbi_tax_db.rds")
