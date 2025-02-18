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
    "ggplot2",
    "httr",
    "kmer",
    "magrittr",
    "methods",
    "openssl",
    "parallel",
    "patchwork",
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
    "seq_names_file",
    "final_db_file",
    "add_root"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# get seq names for key species as list
key_sequences <- readr::read_lines(seq_names_file) %>% as.list()

# species name
spp_name <- key_sequences[[1]] %>% stringr::str_extract(., "[^;]+$")

# import database
final_db <- ape::read.FASTA(final_db_file) %>% DNAbin2DNAstringset()

# emboss ednafull scoring matrix 
ednafull_matrix <- matrix(c(
    # A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
    5, -4, -4, -4, -4,  1,  1, -4, -4,  1, -4, -1, -1, -1, -2,  # A
    -4,  5, -4, -4, -4,  1, -4,  1,  1, -4, -1, -4, -1, -1, -2,  # T
    -4, -4,  5, -4,  1, -4,  1, -4,  1, -4, -1, -1, -4, -1, -2,  # G
    -4, -4, -4,  5,  1, -4, -4,  1, -4,  1, -1, -1, -1, -4, -2,  # C
    -4, -4,  1,  1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1,  # S
    1,  1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1,  # W
    1, -4,  1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1,  # R
    -4,  1, -4,  1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1,  # Y
    -4,  1,  1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1,  # K
    1, -4, -4,  1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1,  # M
    -4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1,  # B
    -1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1,  # V
    -1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1,  # H
    -1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1,  # D
    -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1   # N
), nrow = 15, byrow = TRUE)
# Set row and column names for the EDNAFULL matrix
rownames(ednafull_matrix) <- colnames(ednafull_matrix) <- c("A", "T", "G", "C", "S", "W", "R", "Y", "K", "M", "B", "V", "H", "D", "N")

# parse add_root parameter
if (add_root == "true"){
    add_root <- TRUE
} else {
    add_root <- FALSE
}

if (add_root){
    lineage_ranks <- c("root","kingdom","phylum","class","order","family","genus","species")
} else {
    lineage_ranks <- c("kingdom","phylum","class","order","family","genus","species")
}


### run code

# function that get pairwise percent identity between a single pattern and subject
pairwise_identity <- function(seq_pattern, seq_subject){
    # do alignment
    alignment <- 
        pairwiseAlignment(
            pattern = seq_pattern, 
            subject = seq_subject,
            substitutionMatrix = ednafull_matrix,
            type = "overlap",
            gapOpening = -8, 
            gapExtension = -2,
            scoreOnly = FALSE
        )
  
    # get mismatches against subject
    n_mismatch <- mismatchTable(alignment) %>% nrow() # number of rows in mismatchTable = number of mismatches in overlap
    
    # get number of indel bases against subject
    indels_object <- nindel(alignment)
    n_indels <- insertion(indels_object)[2] + deletion(indels_object)[2]
    
    # combine mismatches and indels 
    n_diff <- n_mismatch + n_indels 
    
    # prop identity of pattern to subject
    # aligned region length: alignedPattern(alignment) %>% lengths
    pattern_identity <- 1 - (n_diff / lengths(alignedPattern(alignment)))
    
    # output identity named with sequence
    names(pattern_identity) <- names(seq_pattern)
    
    return(pattern_identity)
}

## get pid against database for each key sequence
key_pid <- 
    lapply(
        key_sequences,
        function(key, db = final_db, ranks){
        # get pid against database for a single key sequence
        pid_vec <- 
            lapply(
                db,
                pairwise_identity,
                seq_subject = db[names(db) == key]
            ) %>% 
            unlist()
        
        # format tibble of pid and lineage
        pid_tibble <- 
            pid_vec %>% 
            tibble::enframe(value = "pid") %>%
            tidyr::separate(col = name, into = c("id", "lineage_string"), sep = ";", extra = "merge") %>%
            tidyr::separate(col = lineage_string, into = ranks, sep = ";", extra = "merge") %>%
            tidyr::separate(col = id, into = c("seqid", "taxid"), sep = "\\|", extra = "merge") %>%
            dplyr::mutate(
                key = key,
                key_species = stringr::str_extract(key, "[^;]+$")
            )
        
        return(pid_tibble)
        
        },
        db = final_db, 
        ranks = lineage_ranks
    )

# combine pid across key sequences
key_tibble <- 
    key_pid %>%
    dplyr::bind_rows() 

# get non-self sequences identical to key sequences
identical_seqs <- 
    key_tibble %>% 
    dplyr::filter(
        pid == 1, # percent identity is 1
        stringr::str_starts(key, paste0(seqid,"\\|"), negate = TRUE), # seqid is not the same as one of the key sequences
        species != key_species # species is not the same as the key sequence species
    )

# write key tibble to file
readr::write_csv(key_tibble, paste0(spp_name,".key_pid.csv"))

# write identical seqs tibble to file
if ( nrow(identical_seqs) > 0){
    readr::write_csv(identical_seqs, paste0(spp_name, ".identical_seqs.csv"))
} else {
    file.create(paste0(spp_name,".identical_seqs.csv"))
}
