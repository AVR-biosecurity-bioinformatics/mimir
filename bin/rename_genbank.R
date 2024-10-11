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
    "ncbi_rankedlineage_noname"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# get basename chunk for output naming
chunk_index <- tools::file_path_sans_ext(basename(fasta_file)) %>% stringr::str_extract("^.+(?=\\.taxid$)") 

# read in fasta file as DNAbin
seqs <-  ape::read.FASTA(fasta_file)

# read in ncbi tax file
ncbi_rankedlineage_noname <-   readRDS(ncbi_rankedlineage_noname)


### run code

# get names of sequences (.fasta header)
seq_names <- names(seqs)

# split names into accession (seqid) and taxid
seq_names_mat <- stringr::str_split_fixed(seq_names, "\\|", n = 2) 

# name columns
colnames(seq_names_mat) <- c("seqid", "taxid")

# convert matrix to tibble
seq_names_tibble <- 
    tibble::as_tibble(seq_names_mat) %>%
    dplyr::mutate(taxid = as.integer(taxid))

# create valid header format for sequence names
seq_names_header <-
    seq_names_tibble %>%
    # add ncbi taxonomic information per taxid, limited to the allowed ranks
    dplyr::left_join(., ncbi_rankedlineage_noname, by = dplyr::join_by(taxid == tax_id)) %>%
    dplyr::mutate(
        dplyr::across(species:kingdom, .fns = ~replace(., is.na(.), "Unclassified")), # replace NA in ranks columns with "Unclassifed"
        taxid = stringr::str_replace(taxid, "^", "NCBI:") # reformat taxid to have NCBI-specific prefix
    ) %>%
    dplyr::relocate( # reorder columns
        seqid,
        taxid, 
        kingdom,
        phylum,
        class,
        order,
        family,
        genus,
        species
    ) %>%
    tidyr::unite("id", c(seqid, taxid), sep = "|") %>% # combine ids into a single column
    tidyr::unite("ranks", kingdom:species, sep = ";") %>% # combine ranks into a single column
    tidyr::unite("header", id:ranks, sep = ";") # create final header 

# rename sequences with new header
names(seqs) <- seq_names_header$header

# write .fasta to file
write.FASTA(seqs, file = paste0(chunk_index,".renamed.fasta"))
