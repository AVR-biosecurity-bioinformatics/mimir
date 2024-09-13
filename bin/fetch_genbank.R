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
    "taxon",
    "db_file",
    "task_index"
    )
lapply(nf_vars, nf_var_check)

### process variables 

db <- readRDS(db_file)

### run code

## Fetch sequences from GenBank by searching for a taxon name
genbank_seqs <- 
    fetch_seqs(
        x = taxon, 
        database = "genbank", 
        db = db,
        marker="COI[GENE] OR COX1[GENE] OR COXI[GENE]", 
        output = "gb-binom", 
        retry_attempt = 3, 
        retry_wait = 5, 
        multithread = FALSE, 
        quiet=FALSE
    )

# save sequences (DNAbin object) as .rds file
saveRDS(genbank_seqs, paste0(taxon, "_", task_index, "_genbank.rds"))

# write fasta for debugging
if ( params.all_fasta == "true"){
    write_fasta(
        genbank_seqs, 
        file = paste0(taxon, "_", task_index, "_genbank.fasta"), 
        compress = FALSE
        )
}

