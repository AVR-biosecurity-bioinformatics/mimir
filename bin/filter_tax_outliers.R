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
    "cluster_tsv",
    "db_file"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# convert fasta to DNAbin
seqs <- ape::read.FASTA(fasta_file, type = "DNA")

# read cluster TSV 
clusters <- read_delim(cluster_tsv, col_names = c("representative","sequence_name"), delim = "\t")

## params parsing 
cluster_rank <- stringr::str_split_1(params.cluster_rank, ",")
cluster_threshold <- as.numeric(params.cluster_threshold)
cluster_confidence <- as.numeric(params.cluster_confidence)
cluster_nstart <- as.numeric(params.cluster_nstart)
quiet <- FALSE

### run code

# ## find mixed clusters of sequences
# mixed_clusters <- 
#     get_mixed_clusters(
#         x = seqs, # DNAbin list, names must include NCBI taxid
#         db = db, # NCBI taxonomy
#         rank = cluster_rank, # taxonomic rank to check clusters, can be a string or a vector of strings
#         threshold = cluster_threshold, # OTU clustering threshold
#         rngseed = 1, # sets set.seed before clustering to ensure reproducibility
#         return = "consensus", # what to return: "consensus", "all", or "count"
#         k = 5, # k-mer size for input matrix
#         confidence = cluster_confidence, # proportion of sequences that need to have different taxonomy (default 0.8)
#         nstart = cluster_nstart, # random sets chosen for kmeans; higher means more accurate clustering at expense of time 
#         quiet = FALSE
#     ) 

### new-new implementation with MMSeqs2 clustering
# make a tibble of seq names that explicitly encodes order
seq_names <- 
  tibble::tibble(
    sequence_name = names(seqs),
    order = seq(1, length(names(seqs)))
  )
## create named vector from clusters tibble
# get cluster reps
cluster_reps <- clusters$representative %>% unique
# give them new integer names
names(cluster_reps) <- seq(from = 1, to = length(cluster_reps))
# make tibble
cluster_info <- tibble::enframe(cluster_reps, name = "cluster", value = "representative")
rm(cluster_reps)
# join to clusters tibble
clusters_joined <- 
    dplyr::left_join(clusters, cluster_info, by = "representative") %>% 
    dplyr::select(-representative)
rm(clusters)
rm(cluster_info)
# join seq_names to cluster_joined to get otus info
clustered_seqs <- 
    seq_names %>%
    dplyr::left_join(., clusters_joined, by = "sequence_name")
rm(seq_names)
rm(clusters_joined)
# create new named vector ala. old code
otus <- clustered_seqs$cluster
names(otus) <- clustered_seqs$sequence_name
rm(clustered_seqs)
## make sure 'clusters' and 'seqs' are in the same order
# get the taxonomic lineage of each sequence as a tibble/df
lineage <- 
    names(seqs) %>%
    tibble::as_tibble(column_name = "value") %>%
    tidyr::separate(col=value, into=c("id", "lineage_string"), sep=";", extra="merge") %>%
    tidyr::separate(col=lineage_string, into=c("kingdom","phylum", "class", "order", "family", "genus", "species"), sep=";", extra="merge") %>%
    tidyr::separate(col=id, into=c("seqid", "taxid"), sep="\\|", extra="merge") 
# get the taxon names for each sequence at the rank specified
rank_taxa <- lineage %>% pull(!!cluster_rank)
# name the taxon name vector using the sequence ids (ie. accession)
names(rank_taxa) <- lineage$seqid
# remove unneeded
rm(lineage)
# get a vector of clusters + seq names where the rank of interest is not "Unclassified" (unclassified breaks mixed cluster detection)
f <- as.factor(otus[!rank_taxa %>% stringr::str_detect("^Unclassified$")])
rm(otus)
# vector of taxon names at the given rank for classified sequences
rank_classified <- rank_taxa[!rank_taxa %>% stringr::str_detect("^Unclassified$")]
rm(rank_taxa)
# make list where each element is an otu cluster, with the members are the accessions and the taxa 
splitlist <- split(rank_classified, f)
# keep lists where there are more than two members (otherwise can't determine consensus (n = 2) or it is unnecessary (n = 1))
splitlist <- splitlist[tabulate(f) > 2]
# remove unneeded
rm(rank_classified)
rm(f)
# define new function to find mixed clusters, returns a df where "Acc" column is the accession of a potentially bad sequence
find_mixed_new <- function(y) {
    hashes <- paste0(gsub("\\|.*$", "\\1", names(y)), y)
    yu <- y[!duplicated(hashes)]
    if (length(unique(yu)) < 2) {
        return(NULL)
    }
    # Tabulate taxon names
    tab <- sort(table(yu), decreasing = TRUE)
    consensus <- names(tab)[1] # get the consensus taxon name (highest abundance)
    consensus_seqid <- gsub("^.*\\|", "\\1", names(y)[y==consensus][1]) # get a representative member of taxon
    mixed <- y != consensus # potential misannotated (including duplicates)
    mixedu <- yu != consensus # potential misannotated (unique only)
    nu <- length(mixedu) # number of unique sequences
    #Check if there is a clear consensus
    if (tab[1] == tab[2]){
        consensus <- NA
        consensus_seqid <- NA
    }
    res <- 
        data.frame(
            listed = y[mixed],
            consensus = rep(consensus, sum(mixed)),
            consensus_seqid = rep(consensus_seqid,sum(mixed)),
            confidence = sum(!mixedu)/nu,
            cluster_size = nu,
            stringsAsFactors = FALSE
        ) %>%
        tibble::rownames_to_column("Acc")
    return(res)
}
# apply find_mixed_new to list of clusters
mixedtab <- lapply(splitlist, find_mixed_new)
rm(splitlist)
# remove lists where the function returned NULL
mixedtab <- mixedtab[!vapply(mixedtab, is.null, logical(1))]

if ( lapply(mixedtab,nrow) %>% unlist %>% sum == 0 ) { # handle list output or data.frame formats
      if (!quiet) {cat("No mixed clusters at", cluster_rank,   "rank \n")}
    mixed_clusters <- NULL
} else if ( lapply(mixedtab,nrow) %>% unlist %>% sum > 0 ){
    mixedtab <- dplyr::bind_rows(mixedtab, .id="cluster")
    mixedtab <- mixedtab[mixedtab$confidence >= cluster_confidence, ]
    mixedtab <- mixedtab[order(mixedtab$confidence, decreasing = TRUE), ]
    if (!quiet) {cat("identified", length(unique(mixedtab$cluster)), "mixed clusters at", cluster_rank, "rank \n")}
    mixed_clusters <- 
        mixedtab %>%
        dplyr::mutate(
            rank = cluster_rank,
            threshold = cluster_threshold
        ) %>%
        dplyr::mutate_if(is.factor, as.character)
    
    # handle case where mixed_clusters is a df with 0 rows instead of NULL (can't explain code causing this)
    if ( is.data.frame(mixed_clusters) ){
        if (nrow(mixed_clusters) == 0){
            mixed_clusters <- NULL
        }
    }
}

# remove mixed cluster seqs
if (!is.null(mixed_clusters)){
    # Get accession numbers to remove
    rem <- mixed_clusters$Acc
    # create searching pattern
    rem_pattern <- stringr::str_c(rem, collapse = "|")
    # Remove any accessions that were found to be in mixed clusters
    seqs_decontaminated <- seqs[!stringr::str_starts(names(seqs), rem_pattern)]
    # store removed sequenced
    seqs_removed <- seqs[stringr::str_starts(names(seqs), rem_pattern)]
} else {
    seqs_decontaminated <- seqs
    seqs_removed <- NULL
}

# write fasta output
write_fasta(
    seqs_decontaminated, 
    file = "seqs_decontaminated.fasta", 
    compress = FALSE
)

if ( is.null(seqs_removed) ){
    file.create("removed.fasta")
} else {
    write_fasta(
        seqs_removed, 
        "removed.fasta"
    )
}
