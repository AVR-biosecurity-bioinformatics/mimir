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
    "primer_fwd",
    "primer_rev"
    )
lapply(nf_vars, nf_var_check)

### process variables 

### run code

# get sequences of all four primers
fwd_ag.string <- primer_fwd
rev_ag.string <- primer_rev
fwd_ag.seq <- Biostrings::DNAString(fwd_ag.string)
rev_ag.seq <- Biostrings::DNAString(rev_ag.string)
fwd_rc.seq <- Biostrings::reverseComplement(fwd_ag.seq)
rev_rc.seq <- Biostrings::reverseComplement(rev_ag.seq)

## disambiguated sequences, named with degenerate primer name and a number from 1 to length()
# fwd_ag
fwd_ag.da <- fwd_ag.seq %>% DNAStringSet(.) %>% DECIPHER::Disambiguate(.) %>% .[[1]]
names(fwd_ag.da) <- seq(1,length(fwd_ag.da)) %>% stringr::str_replace(., "^", "fwd_ag")

# rev_ag
rev_ag.da <- rev_ag.seq %>% DNAStringSet(.) %>% DECIPHER::Disambiguate(.) %>% .[[1]]
names(rev_ag.da) <- seq(1,length(rev_ag.da)) %>% stringr::str_replace(., "^", "rev_ag")

# fwd_ag
fwd_rc.da <- fwd_rc.seq %>% DNAStringSet(.) %>% DECIPHER::Disambiguate(.) %>% .[[1]]
names(fwd_rc.da) <- seq(1,length(fwd_rc.da)) %>% stringr::str_replace(., "^", "fwd_rc")

# rev_rc
rev_rc.da <- rev_rc.seq %>% DNAStringSet(.) %>% DECIPHER::Disambiguate(.) %>% .[[1]]
names(rev_rc.da) <- seq(1,length(rev_rc.da)) %>% stringr::str_replace(., "^", "rev_rc")

# write to .fasta
c(fwd_ag.da, rev_ag.da, fwd_rc.da, rev_rc.da) %>% ape::as.DNAbin() %>% ape::write.FASTA(., "primers.fasta")