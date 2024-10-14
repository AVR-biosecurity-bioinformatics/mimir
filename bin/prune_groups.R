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
    "seqs_file"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read sequences from file
seqs <- readRDS(seqs_file)

## parse params
max_group_size <- as.numeric(params.max_group_size)

prune_method <- params.prune_method

if (params.remove_unclassified == "true") {
    remove_unclassified <- TRUE
} else {
    remove_unclassified <- FALSE
} 

### run code

## prune groups of sequences with identical taxonomic IDs down to a certain number

# ### old function
# seqs_pruned <- 
#     prune_groups(
#         x = seqs, # DNAbin or DNAStringset object
#         max_group_size = max_group_size, # max sequences to keep
#         discardby = prune_method, # 'length': discard smallest sequences first; 'random': discard randomly
#         dedup = TRUE, # remove sequences with identical taxonomic name and sequence first
#         prefer = NULL, # vector of sequence names to prefer (eg. high-quality internal sequences)
#         quiet = FALSE
#     )

### new function
prune_groups_alt <- function(x, max_group_size = 5, dedup = TRUE, discardby = "length", prefer=NULL, quiet = FALSE, remove_unclassified = FALSE) {
  # Convert to DNAbin
  if (!methods::is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
    if (all(is.na(ape::base.freq(x)))) {stop("Error: Object is not coercible to DNAbin \n")}
  }
  if (!is.null(prefer) & !is.character(prefer)){
    stop("prefer must be either NULL or a vector of sequence names to prefer")
  }

  # Remove duplicate sequences
  if (dedup) {
    dup <- length(x)
    # gets species-level ID only
    taxids <- names(x) %>% stringr::str_remove("^.*;")

    ## Consider taxonomic name and sequence identity in deduplication
    # get hash for each name and sequence together
    hashes <- purrr::map2(x, taxids, ~{
      openssl::md5(paste(c(.y, as.vector(.x)), collapse=""))
    }) %>%
      unlist()

    # get list of all duplicated hashes
    dupes <- unique(hashes[duplicated(hashes)]) # vector of duplicated hashes
    remove <- logical(length(x))
    for (i in 1:length(dupes)) {
      index <- which(hashes %in% dupes[i])
      # Check if in prefered - If so keep first element in the prefered
      if (is.character(prefer) & any(names(x)[index] %in% prefer)){
        keep <- index[names(x)[index] %in% prefer][1]
      } else {
        # otherwise just pick first element
        keep <- index[1]
      }
      remove[index[!index %in% keep]] <- TRUE
    }
    x <- x[!remove]
    if (!quiet) cat(paste0((dup - length(x)), " duplicate sequences removed \n"))
  }

  ## remove sequences with "Unclassified" within its lineage string
  #### TODO: make sure this checks for preferred sequences first and doesn't remove them
  if (remove_unclassified) {
    # get lineage string
    lineage_string <- names(x) %>% stringr::str_remove("^.+?;")
    # find lineage strings that contain "Unclassified"
    remove <- lineage_string %>% stringr::str_detect("Unclassified")
    # remove sequences
    x <- x[!remove]
    if (!quiet) cat(paste0(sum(remove), " sequences pruned with unclassified taxonomy \n"))

  }
  
  # Remove sequences from groups where more species names than max_group_size
  # "acc" is considered seqid + taxid, while "taxon" is the lineage string
  groups <- names(x) %>%
    stringr::str_split_fixed(";", n = 2) %>%
    tibble::as_tibble(.name_repair = ~ c("acc", "taxon")) %>%
    dplyr::pull(taxon) %>%
    openssl::md5()
  groupCounts <- table(groups) # Count number of seqs per group
  u_groups <- names(groupCounts) # Get unique groups

  remove <- logical(length(x))
  # Random discarding
  if (discardby == "random") {
    for (i in which(groupCounts > max_group_size)) {
      index <- which(groups == u_groups[i])
      # Deal with prefered
      if (is.character(prefer) & any(names(x)[index] %in% prefer)){
        in_pref <- index[names(x)[index] %in% prefer]
        # Check if can sample just from prefered, or need a mix of the two
        if(in_pref >= max_group_size){
          keep <- sample(
            x = in_pref,
            size = max_group_size,
            replace = FALSE
          )
        } else{
          to_sample <- max_group_size - length(in_pref)
          keep <- c(in_pref,
                    sample(
                      x = index[-in_pref],
                      size = to_sample,
                      replace = FALSE
                    ))
        }
      } else {
        # otherwise just take random sample
        keep <- sample(
          x = index,
          size = max_group_size,
          replace = FALSE
        )
      }
      remove[index[-keep]] <- TRUE
    }
    # length discarding
  } else if (discardby == "length") {
    for (i in which(groupCounts > max_group_size)) {
      index <- which(groups == u_groups[i])
      rem <- sapply(x[index], function(s) length(s) - sum(s == as.raw(c(4)))) # Get lengths
      names(rem) <- index
      rem <- sort(rem, decreasing = TRUE)

      # deal with prefered
      if (is.character(prefer) & any(names(x)[index] %in% prefer)){
        in_pref <- index[names(x)[index] %in% prefer]
        non_pref <- index[!names(x)[index] %in% prefer]

        # Get the minimum length of sequences that are being dropped
        min_len <- min(rem[1:max_group_size])

        # Check if there are mixed prefered and non prefered at that length
        mixed_pref <- rem[names(rem) %in% in_pref] == min_len
        mixed_pref <- names(mixed_pref)[mixed_pref]
        mixed_non_pref <- rem[names(rem) %in% non_pref] == min_len
        mixed_non_pref <- names(mixed_non_pref)[mixed_non_pref]

        if((length(mixed_pref) > 0) & (length(mixed_pref) > 0)){
          # re-sort vector to prefer
          sample_from <- as.integer(c(names(rem)[rem > min_len], mixed_pref, mixed_non_pref, names(rem)[rem < min_len]))
          keep <- as.integer(sample_from[1:max_group_size])
        } else {
          # otherwise just take the longest
          keep <- as.integer(names(rem[1:max_group_size]))
        }
      } else{
        keep <- as.integer(names(rem[1:max_group_size]))
      }
      remove[index[!index %in% keep]] <- TRUE
    }
  }
  x <- x[!remove]
  if (!quiet) cat(paste0(sum(remove), " sequences pruned from over-represented groups \n"))
  return(x)
}

# using new function
seqs_pruned <- 
    prune_groups_alt(
        x = seqs, # DNAbin or DNAStringset object
        max_group_size = max_group_size, # max sequences to keep
        discardby = prune_method, # 'length': discard smallest sequences first; 'random': discard randomly
        dedup = TRUE, # remove sequences with identical taxonomic name and sequence first
        prefer = NULL, # vector of sequence names to prefer (eg. high-quality internal sequences)
        quiet = FALSE,
        remove_unclassified = remove_unclassified # remove sequences with any unclassified taxonomic rank
    )


# save filtered sequences as .rds file
saveRDS(seqs_pruned, "seqs_pruned.rds")

# write fasta
write_fasta(
    seqs_pruned, 
    file = "seqs_pruned.fasta", 
    compress = FALSE
)
