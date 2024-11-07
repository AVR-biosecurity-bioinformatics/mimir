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
    "fasta_file"
)
lapply(nf_vars, nf_var_check)

### process variables 

# read sequences from file
seqs <- ape::read.FASTA(fasta_file)

# params
orient <- FALSE
max_group_size <- as.numeric(params.idtaxa_max_group_size)
max_iterations <- as.numeric(params.idtaxa_max_iterations)
if ( params.idtaxa_allow_group_removal == "true" ){
    allow_group_removal <- TRUE
} else {
    allow_group_removal <- FALSE
}
taxid <- NULL # null because we have lineage in name
quiet <- FALSE

# Convert to DNAstringset for DECIPHER
if(methods::is(seqs, "DNAbin")){
  seqs <-  DNAbin2DNAstringset(seqs)
}

# Remove gaps
if(!quiet){ message("Removing gaps")}
seqs <- DECIPHER::RemoveGaps(seqs)

# orient if requested
if(orient){
  if(!quiet){ message("Orienting sequences")}
  seqs <- DECIPHER::OrientNucleotides(seqs)
}


#Add Root Rank if not already there (TODO: make this code check if present directly rather than relying on parameter value)
if ( params.add_root == "false" ){
    names(seqs) <- names(seqs)  %>% stringr::str_replace(";[;]*", ";Root;")
}

# obtain the taxonomic assignments
groups <- names(seqs) # sequence names
# assume the taxonomy begins with 'Root;'
groups <- gsub("(.*)(Root;)", "\\2", groups) # extract the group label
groupCounts <- table(groups)
u_groups <- names(groupCounts) # unique groups

if (!quiet) {message(paste0(length(u_groups), " unique species"))}

# Pruning training set
remove <- logical(length(seqs))
for (i in which(groupCounts > max_group_size)) {
  index <- which(groups==u_groups[i])
  keep <- sample(length(index),
                 max_group_size)
  remove[index[-keep]] <- TRUE
}
if (!quiet) {message(paste0(sum(remove), " sequences pruned from over-represented groups"))}

# Training the classifier
probSeqsPrev <- integer() # suspected problem sequences from prior iteration

train <- seqs[!remove]
taxonomy <- gsub("(.*)(Root;)", "\\2", names(train))
remove <- logical(length(train))

for (i in seq_len(max_iterations)) {
  if (!quiet) {message("Training iteration: ", i, "\n", sep="")}
  # train the classifier
  trainingSet <- DECIPHER::LearnTaxa(train[!remove], taxonomy[!remove], taxid)

  # look for problem sequences
  probSeqs <- trainingSet$problemSequences$Index
  if (length(probSeqs)==0) {
    if (!quiet) {message("No problem sequences remaining.\n")}
    break
  } else if (length(probSeqs)==length(probSeqsPrev) &&
             all(probSeqsPrev==probSeqs)) {
    if (!quiet) {message("Iterations converged.\n")}
    break
  }
  if (i==max_iterations)
    break
  probSeqsPrev <- probSeqs

  # remove any problem sequences
  index <- which(!remove)[probSeqs]
  remove[index] <- TRUE # remove all problem sequences
  if (!allow_group_removal) {
    # replace any removed groups
    missing <- !(u_groups %in% groups[!remove])
    missing <- u_groups[missing]
    if (length(missing) > 0) {
      index <- index[groups[index] %in% missing]
      remove[index] <- FALSE # don't remove
    }
  }
}

if (!quiet) {message(paste0(sum(remove), " sequences removed"))}
if (!quiet) {message(paste0(length(probSeqs), " problem sequences remaining"))}

# save model
saveRDS(trainingSet, "idtaxa_model.rds")