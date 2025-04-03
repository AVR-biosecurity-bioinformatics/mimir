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
    "gb_file",
    "accessions_file",
    "ncbi_rankedlineage_noname",
    "placeholder_as_unclassified"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read in sequences GenBank flat file
gb_raw <- scan(file = gb_file, what = "", sep = "\n", quiet = TRUE)

# read in accessions file
accessions <- scan(file = accessions_file, what = "", sep = "\n", quiet = TRUE)

# read in ncbi tax file
ncbi_rankedlineage_noname <-   readRDS(ncbi_rankedlineage_noname)

# parse params.placeholder_as_unclassified
if ( placeholder_as_unclassified == "true" ){
    placeholder_as_unclassified <- TRUE
} else if ( placeholder_as_unclassified == "false" ){
    placeholder_as_unclassified <- FALSE
} else {
    stop("'placeholder_as_unclassified' is not 'true' or 'false'")
}

### run code

# Truncate record at last // to avoid broken records
gb <- tryCatch(
    { gb_raw[1:max(which(grepl("^//", gb_raw)))] }, 
    error = function(e){
        warning("Failed parsing")
        print(max(which(grepl("^//", gb_raw))))
        return(NULL)
    }
)

# check for accessions
if(sum(grepl("ACCESSION", gb)) < 1){
    stop("No accessions in GenBank flat file")
}

# lines where records start and stop
record_start <- which(grepl("^LOCUS +", gb))
record_stop <- which(grepl("^//", gb))

# check same number of start and stop positions
if (length(record_start) == length(record_stop)){
    n_seqs <- length(record_start)
} else {
    stop("Incorrect length of start and stop")
}

# extract records into a list
records <- vector("list", length = n_seqs)
for (l in 1:n_seqs){
    records[[l]] <- gb[(record_start[l]):(record_stop[l])]
}

# function to format GB record into named sequence
sequence_from_gb_record <- function(x){
    # get versions (accession with version decimal)
    seq_ver <- gsub("+VERSION +", "", grep("VERSION", x, value = TRUE))
    # get taxids
    seq_taxid <- grep("db_xref=\\\"taxon:", x, value = TRUE) %>% stringr::str_remove_all("( +/db_xref=\"taxon:)|(\\\")")
    # get "ORIGIN" line (just before nucleotide sequence)
    seq_start <- which(grepl("^ORIGIN", x))
    # get "//" line just after sequence
    seq_end <- which(grepl("^//", x))
    # check if there is a sequence
    if ( length(seq_start) == 0 || length(seq_end) == 0 ){
        print(paste0("No sequence available for ", seq_ver))
        return(NULL)
    } else {
        # get sequence
        sequence <- toupper(paste(stringr::str_remove_all(x[(seq_start+1):(seq_end-1)], "[^A-Za-z]"), collapse=""))
    }
    # check only one version, taxid and sequence
    if(length(seq_ver) > 1){
        stop(paste0("More than one sequence version:", stringr::str_flatten_comma(seq_ver)))
    }
    if(length(seq_taxid) > 1){
        if(seq_taxid %>% unique %>% length == 1){ # if taxid is mentioned more than once...
            seq_taxid <- unique(seq_taxid) # collapse if all the same
        } else {
            print(x)
            stop(paste0("More than one unique sequence taxid:", stringr::str_flatten_comma(seq_taxid),"for seq_ver ",seq_ver))
        }
    }
    if(length(sequence) > 1){
        stop(paste0("More than one sequence string:", stringr::str_flatten_comma(sequence)))
    }
    # check sequence contains only valid nucleotide characters
    if ( sequence %>% stringr::str_detect("[^GATCRYMKSWHBVDN]") ){
        stop(paste0("Sequence string for ", seq_ver, " contains non-IUPAC characters"))
    }
    # name sequence
    names(sequence) <- paste0(seq_ver,"|",seq_taxid)
    print(paste0("Parsed sequence ",seq_ver))
    return(sequence)
}

# for each element of list (record), get version, taxid and sequence and put into DNAbin object
seqs <- 
    lapply(
        records, 
        sequence_from_gb_record
    ) %>%
    # remove null elements (no sequence)
    purrr::discard(is.null) %>%
    unlist() %>%
    char2DNAbin()

# write list of accessions not parsed into sequences
seq_accessions <- names(seqs) %>% stringr::str_extract(., "^.+(?=\\|)")
noseq_accessions <- accessions[!accessions %in% seq_accessions] 
if ( length(noseq_accessions) > 0 ){
    write_lines(noseq_accessions, "accessions_failed.txt")
} else {
    file.create("accessions_failed.txt")
}

## add taxonomic lineage string
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
    # conditionally replace 'placeholder' species names (eg. "Genus sp. XYZ") with "Unclassified"
    {
        if (placeholder_as_unclassified) {
            dplyr::mutate(
                .,
                species = dplyr::case_when(
                    stringr::str_detect(species, paste0(genus, " [:alnum:]+\\.( |$)")) ~ "Unclassified",
                    .default = species 
                )
            )
        } else { . }
    } %>%
    dplyr::mutate(
        dplyr::across(species:kingdom, .fns = ~replace(., is.na(.), "Unclassified")), # replace NA in ranks columns with "Unclassifed"
        # dplyr::across(species:kingdom, .fns = ~stringr::str_replace_all(., "[ \\/:\\(\\)&,'#<>]", "_")), # replace problematic characters in lineage string with underscores
        # dplyr::across(species:kingdom, .fns = ~stringr::str_replace_all(., "_+", "_")), # replace two or more underscores in a row with a single underscore in lineage string
        # dplyr::across(tidyselect::everything(), .fns = ~stringr::str_replace_all(., " ", "_")), # replace all spaces with underscores
        dplyr::across(species:kingdom, .fns = ~stringr::str_replace_all(ranks, " +", " ")), # replace runs of spaces with a single space (to allow HMMER output parsing later)
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
write.FASTA(seqs, file = "renamed.fasta")
