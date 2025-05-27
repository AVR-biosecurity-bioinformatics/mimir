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
    "fasta_files",
    "counts_file",
    "dist_threshold"
    )
lapply(nf_vars, nf_var_check)

### process variables 

## read in list of sequences
fasta_vector <- # convert Groovy to R vector format
    stringr::str_extract_all(fasta_files, pattern = "[^\\s,\\[\\]]+") %>% unlist()

# read into list
seqs_list <- lapply(fasta_vector, ape::read.FASTA) 

# name DNAbin objects in list based on .fasta names
names(seqs_list) <- fasta_vector %>% stringr::str_remove(., "\\.retained\\.aligned\\.fasta") %>% stringr::str_remove(., "^.*/")

# check all files are .fasta 
for (i in 1:length(fasta_vector)) {
    if(!stringr::str_detect(fasta_vector[i], "\\.fa(sta)$")){stop(paste0("Sequence file '",fasta_vector[i],"' does not appear to be a FASTA file"))}
}

## read counts .csv into tibble
counts_tibble <- 
    readr::read_tsv(counts_file, col_names = c("name", "count")) %>%
    dplyr::mutate(
        name = stringr::str_remove(name, "^>")
    )

## parse params
dist_threshold <- as.numeric(dist_threshold)

### run code

## function to detect sequence outliers within a species
remove_intraspp_dist_outliers <- function(seqs, removal_threshold = 0.05, counts) {
    # convert sequences
    seqs <- seqs %>% DNAbin2DNAstringset(remove_gaps = FALSE) 
    # get species name
    spp_name <-  
        seqs %>%
        names() %>%
        tibble::as_tibble(column_name = "value") %>%
        tidyr::separate(col=value, into=c("id", "lineage_string"), sep=";", extra="merge") %>%
        tidyr::separate(col=lineage_string, into = c("kingdom","phylum", "class", "order", "family", "genus", "species"), sep=";", extra="merge") %>%
        tidyr::separate(col=id, into=c("seqid", "taxid"), sep="\\|", extra="merge") %>%
        dplyr::pull(species) %>%
        unique()
    # check only one species name
    if (length(spp_name) > 1) {stop(paste0("Multiple species detected when one needed: '",stringr::str_c(spp_name, sep = ","),"'"))}
    # calculate intraspecific distance if...
    if (
        spp_name != "Unclassified" && # species is not Unclassified
        length(unique(seqs)) > 1  # there are at least 2 unique sequences
    ) {
        # check sequences are aligned by comparing lengths
        if(length(table(lengths(seqs))) > 1){ stop(paste0("Sequences for species '",spp_name,"' don't appear to be aligned"))} 
        # check seq names are in counts tibble
        if ( !all(names(seqs) %in% counts_tibble$name) ){
            stop(paste0("Not all sequence names for '",spp_name,"' are in counts_tibble"))
        }
        # get sequences as a tibble
        seqs_tibble <- 
            seqs %>% 
            as.data.frame %>% 
            tibble::as_tibble(rownames = "name") %>%
            dplyr::rename(seq = x) %>%
            # join seqs to counts from FILTER_REDUNDANT
            dplyr::left_join(., counts_tibble, by = "name")
        # get vector of sequence weights (counts)
        seqs_weights <- seqs_tibble$count
        # convert back to DNAStringset
        seqs_reps <- tibble::deframe(seqs_tibble[,1:2]) %>% DNAStringSet(.)
        # get distance matrix for representative sequences
        message(paste0("Calculating distance matrix for species '",spp_name,"'"))
        distmat <- 
          DECIPHER::DistanceMatrix(
            seqs_reps, 
            includeTerminalGaps = FALSE, 
            penalizeGapLetterMatches = TRUE,
            #penalizeGapGapMatches = FALSE, 
            correction = "Jukes-Cantor",
            processors = 1, 
            verbose = TRUE
          )
        # find index of most central element, using weighted median
        # code from: https://stackoverflow.com/questions/2748725/is-there-a-weighted-median-function
        # works because weights are counts and therefore integers
        center <- which.min(
            apply(
                distmat,
                1,
                function(x, w){
                    median(rep(x, times = w))
                },
                w = seqs_weights
            )
        )
        gc() # clean memory
        # distance from each element to most central element
        dd <- distmat[center,]
        # remove distance matrix from memory
        rm(distmat)
        gc()
        # give distance to sequences
        seqs_dd <- 
            seqs_tibble %>%
            dplyr::mutate(dd = dd) 
        if (!all(seqs_dd$name == names(seqs))){stop("Seqs not in the correct order when subsetting by distance")}
        # remove sequences above threshold
        if(any(seqs_dd$dd > removal_threshold)){
            removed <- ape::as.DNAbin(seqs[seqs_dd$dd > removal_threshold]) %>% ape::del.gaps() # remove alignment
            retained <- ape::as.DNAbin(seqs[seqs_dd$dd <= removal_threshold]) %>% ape::del.gaps() # remove alignment
        } else {
            removed <- NULL
            retained <- ape::as.DNAbin(seqs) %>% ape::del.gaps() # remove alignment
        }
    } else {
        # if sequences are not to be compared, save them all as retained
        removed <- NULL
        retained <- ape::as.DNAbin(seqs) %>% ape::del.gaps() # remove alignment
    }
    # output list of removed and retained sequences
    out <- list(removed, retained)
    names(out) <- c("removed","retained")
    message(paste0("Finished species '",spp_name,"'"))
    gc() # clean memory
    return(out)
}

## loop through each DNAbin in seqs_list, calculating distance matrix and saving removed and retained sequences as .fasta files
purrr::imap(
    .x = seqs_list,
    .f = \(
        x, # element of seqs_list
        idx, # name of element (ie. lineage name)
        removal_threshold = dist_threshold
    ){
    # get distance matrix and output removed and retained sequences
    intraspp_out <- remove_intraspp_dist_outliers(x, removal_threshold, counts = counts_tibble)
    gc() # clean memory if needed
    # write sequences to file
    if (!is.null(intraspp_out$removed)){
        suppressMessages(
            write_fasta(
                intraspp_out$removed,
                file = paste0(idx,".removed.fasta")
            )
        )
    } else {
      file.create(paste0(idx,".removed.fasta"))
    }
    # save retained seq as .fasta (empty file if not present)
    if (!is.null(intraspp_out$retained)){
      suppressMessages(
            write_fasta(
                intraspp_out$retained,
                file = paste0(idx,".retained.fasta")
            )
      )
    } else {
      file.create(paste0(idx,".retained.fasta"))
    }
    message(paste0("Removed (",intraspp_out$removed %>% length,") and retained (",intraspp_out$retained %>% length,") sequences from '",idx,"' saved to file"))
  }
)


