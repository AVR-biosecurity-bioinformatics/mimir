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
    "db_file"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# convert fasta to DNAbin
seqs <- ape::read.FASTA(fasta_file, type = "DNA")

# read ncbi db file
db <- readRDS(db_file)

### run code

## find mixed clusters of sequences
mixed_clusters <- 
    get_mixed_clusters(
        x = seqs, # DNAbin list, names must include NCBI taxid
        db = db, # NCBI taxonomy
        rank = "species", # taxonomic rank to check clusters, can be a string or a vector of strings
        threshold = 0.97, # OTU clustering threshold
        rngseed = 1, # sets set.seed before clustering to ensure reproducibility
        return = "consensus", # what to return: "consensus", "all", or "count"
        k = 5, # k-mer size for input matrix
        confidence = 0.6, # proportion of sequences that need to have different taxonomy (default 0.8)
        nstart = 20, # random sets chosen for kmeans; higher means more accurate clustering at expense of time 
        quiet = FALSE
    ) 

# Get accession numbers to remove
rem <- mixed_clusters$Acc

# Get current accession numbers
acc <- str_replace(names(seqs), "(?:.(?!;))+$", "")

# Remove any accessions that were found to be in mixed clusters
seqs_decontaminated <- seqs[!acc %in% rem]

# save filtered sequences as .rds file
saveRDS(seqs_decontaminated, "seqs_decontaminated.rds")

# write fasta for debugging
if ( params.all_fasta == "true"){
    write_fasta(
        seqs_decontaminated, 
        file = "seqs_decontaminated.fasta", 
        compress = FALSE
        )
}