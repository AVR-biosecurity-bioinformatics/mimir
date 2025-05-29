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
    "coding",
    "marker_type",
    "ncbi_gencodes"
)
lapply(nf_vars, nf_var_check)

### process variables 

# read in seqs from file
seqs <- ape::read.FASTA(fasta_file)

# import ncbi_gencodes
ncbi_gencodes <- readRDS(ncbi_gencodes)

# check validity of marker_type
if (!marker_type %in% c("nuclear","mitochondrial","plastid","hydrogen")){
    stop(paste0("ERROR: 'marker_type' is '", marker_type,"' but must be one of 'nuclear', 'mitochondrial', 'plastid' or 'hydrogen'"))
}

if (coding == "true"){
    coding <- TRUE
} else if (coding == "false") {
    coding <- FALSE
} else {
    stop("Invalid 'coding' variable")
}

allowed_ranks <- c("kingdom","phylum", "class", "order", "family", "genus", "species")

### run code

# deduplicate sequences
seqs_dss <- 
    seqs %>% 
    DNAbin2DNAstringset() %>%
    as.character(use.names = TRUE) %>%
    tibble::enframe() %>%
    dplyr::group_by(name, value) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    tibble::deframe() %>%
    Biostrings::DNAStringSet()

# get the taxonomic lineage of each sequence as a tibble/df
lineage <- 
    names(seqs_dss) %>%
    tibble::as_tibble(column_name = "value") %>%
    tidyr::separate(col=value, into=c("id", "lineage_string"), sep=";", extra="merge") %>%
    tidyr::separate(col=lineage_string, into=allowed_ranks, sep=";", extra="merge") %>%
    tidyr::separate(col=id, into=c("seqid", "taxid"), sep="\\|", extra="merge") 

## match gencodes based on lowest NCBI rank, progressively censoring ranks until a match is found
# match gencodes based on all ranks
match_all <- 
    lineage %>%
    dplyr::mutate(row_order = row_number()) %>%  # make a new column that contains original row order
    dplyr::left_join(., ncbi_gencodes, by = dplyr::join_by(kingdom, phylum, class, order, family, genus, species)) %>% 
    dplyr::filter(!is.na(tax_id))

# censor species and match
match_censor_s <- 
    lineage %>%
    dplyr::mutate(
        row_order = row_number(),
        species = "Unclassified"
        ) %>%
    dplyr::left_join(., ncbi_gencodes, by = dplyr::join_by(kingdom, phylum, class, order, family, genus, species)) %>% 
    dplyr::filter(
        !is.na(tax_id), 
        !seqid %in% match_all$seqid
        )

# censor genus and match
match_censor_g <- 
    lineage %>%
    dplyr::mutate(
        row_order = row_number(),
        species = "Unclassified",
        genus = "Unclassified"
    ) %>%  # make a new column that contains original row order
    dplyr::left_join(., ncbi_gencodes, by = dplyr::join_by(kingdom, phylum, class, order, family, genus, species)) %>% 
    dplyr::filter(
        !is.na(tax_id), 
        !seqid %in% c(match_all$seqid,match_censor_s$seqid)
    )

# censor family and match
match_censor_f <- 
    lineage %>%
    dplyr::mutate(
        row_order = row_number(),
        species = "Unclassified",
        genus = "Unclassified",
        family = "Unclassified"
    ) %>%  # make a new column that contains original row order
    dplyr::left_join(., ncbi_gencodes, by = dplyr::join_by(kingdom, phylum, class, order, family, genus, species)) %>% 
    dplyr::filter(
        !is.na(tax_id), 
        !seqid %in% c(match_all$seqid,match_censor_s$seqid, match_censor_g$seqid)
    )

# censor order and match
match_censor_o <- 
    lineage %>%
    dplyr::mutate(
        row_order = row_number(),
        species = "Unclassified",
        genus = "Unclassified",
        family = "Unclassified",
        order = "Unclassified"
    ) %>%  # make a new column that contains original row order
    dplyr::left_join(., ncbi_gencodes, by = dplyr::join_by(kingdom, phylum, class, order, family, genus, species)) %>% 
    dplyr::filter(
        !is.na(tax_id), 
        !seqid %in% c(match_all$seqid,match_censor_s$seqid, match_censor_g$seqid, match_censor_f$seqid)
    )


# censor class and match
match_censor_c <- 
    lineage %>%
    dplyr::mutate(
        row_order = row_number(),
        species = "Unclassified",
        genus = "Unclassified",
        family = "Unclassified",
        order = "Unclassified",
        class = "Unclassified"
    ) %>%  # make a new column that contains original row order
    dplyr::left_join(., ncbi_gencodes, by = dplyr::join_by(kingdom, phylum, class, order, family, genus, species)) %>% 
    dplyr::filter(
        !is.na(tax_id), 
        !seqid %in% c(match_all$seqid,match_censor_s$seqid, match_censor_g$seqid, match_censor_f$seqid, match_censor_o$seqid)
    )

# censor phylum and match
match_censor_p <- 
    lineage %>%
    dplyr::mutate(
        row_order = row_number(),
        species = "Unclassified",
        genus = "Unclassified",
        family = "Unclassified",
        order = "Unclassified",
        class = "Unclassified",
        phylum = "Unclassified"
    ) %>%  # make a new column that contains original row order
    dplyr::left_join(., ncbi_gencodes, by = dplyr::join_by(kingdom, phylum, class, order, family, genus, species)) %>% 
    dplyr::filter(
        !is.na(tax_id), 
        !seqid %in% c(match_all$seqid,match_censor_s$seqid, match_censor_g$seqid, match_censor_f$seqid, match_censor_o$seqid, match_censor_c$seqid)
    )

# censor kingdom and match
match_censor_k <- 
    lineage %>%
    dplyr::mutate(
        row_order = row_number(),
        species = "Unclassified",
        genus = "Unclassified",
        family = "Unclassified",
        order = "Unclassified",
        class = "Unclassified",
        phylum = "Unclassified",
        kingdom = "Unclassified"
    ) %>%  # make a new column that contains original row order
    dplyr::left_join(., ncbi_gencodes, by = dplyr::join_by(kingdom, phylum, class, order, family, genus, species)) %>% 
    dplyr::filter(
        !is.na(tax_id), 
        !seqid %in% c(match_all$seqid,match_censor_s$seqid, match_censor_g$seqid, match_censor_f$seqid, match_censor_o$seqid, match_censor_c$seqid, match_censor_p$seqid)
    )

# combine matched rows together and arrange by original order
matched_gencode <- 
    dplyr::bind_rows(
        match_all,
        match_censor_s,
        match_censor_g,
        match_censor_f,
        match_censor_o,
        match_censor_c,
        match_censor_p,
        match_censor_k
    ) %>%
    dplyr::select(seqid, taxid, row_order, gencode, mtgencode, plgencode, hdgencode) %>%
    dplyr::arrange(row_order) # keep in same order as input sequences

# get vector of genetic codes for each sequence depending on marker type
genetic_code_v <- 
    switch(
        marker_type,
        nuclear = matched_gencode$gencode,
        mitochondrial = matched_gencode$mtgencode,
        plastid = matched_gencode$plgencode,
        hydrogen = matched_gencode$hdgencode
    )

# check length of seqs is the same as genetic_code_v
if ( length(seqs_dss) != length(genetic_code_v) ){
    stop("Number of sequences in input .fasta does not match number of returned genetic codes")
}

# clean up memory
rm(lineage)
rm(match_all)
rm(match_censor_s)
rm(match_censor_g)
rm(match_censor_f)
rm(match_censor_o)
rm(match_censor_c)
rm(match_censor_p)
rm(match_censor_k)
gc()

# create frame for each sequence (one list per frame)
frames <- 
    c(
        lapply(1:3, function(pos) XVector::subseq(seqs_dss, start=pos)),
        lapply(1:3, function(pos) XVector::subseq(Biostrings::reverseComplement(seqs_dss), start=pos))
    ) 
 
names(frames) <- c("1","2","3",".1",".2",".3")

# convert list of frames (6 elements) to list of sequence frames (n sequences elements)
frames_list <- 
    frames %>%
    purrr::imap(
        .x = ., 
        function(x, idx){
            tibble::tibble(
                name = names(x),
                sequence = as.character(x),
                frame = idx
            ) %>%
            dplyr::mutate(
                frame = stringr::str_replace(frame, "\\.","-") # replace '.' with '-'
            )
        }
    ) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct() %>% # remove possibly exactly duplicated sequences
    split(., .$name) %>%
    # convert back to DNAStringset format, per input sequence
    lapply(
        ., 
        function(x) {
            x_df <- as.data.frame(x)
            seqs_DSS <- DNAStringSet(x_df$sequence)
            names(seqs_DSS) <- x$name %>% stringr::str_replace(., "$", paste0("_frame=",x$frame)) # add frame to name
            return(seqs_DSS)
        } 
    )

# check identical length to genetic code vector
if ( length(frames_list) != length(genetic_code_v) ) {
    stop("'frames_list' and 'genetic_code_v' are of different lengths when they should be equal!")
}

# translate frames using genetic codes vector
translated <- 
    mapply( # applies function using combination of two vectors 
        function(y, z) {
            suppressWarnings(
                Biostrings::translate(
                    x = y, 
                    genetic.code = Biostrings::getGeneticCode(z),
                    if.fuzzy.codon=c("solve", "X"),
                    no.init.codon = TRUE # don't treat each sequence as starting with a start codon
                )
            )
        }, 
        frames_list, # y
        genetic_code_v # z
    ) %>%
    # get positions of stop codons within each sequence and append to name 
    lapply(
        ., 
        function(x){
            sc_positions <- 
                x %>%
                # get amino acid residues as character string
                as.character() %>% 
                # get locations of all stop codons as a matrix
                stringr::str_locate_all(., "\\*") %>% 
                # select only the first column of matrix, convert to numeric vector, then collapse to string separated by "|"
                lapply(., function(x){
                    as.integer(x[,1]) %>% 
                        stringr::str_flatten(collapse = "|") %>% 
                        stringr::str_replace(., "^", "_sc=")
                }
                ) %>% 
                unlist()
            
            names(x) <- names(x) %>% stringr::str_replace(., "$", sc_positions)
            
            return(x)
        }
    )

rm(frames_list)
gc()

# create a single object from the list of translated sequences
translations_combined <- AAStringSet()
for (i in 1:length(translated)){
  translations_combined <- S4Vectors::append(translations_combined, translated[[i]], after = length(translations_combined))
}

rm(translated)
gc()

# convert stop codons to X character
translations_combined <- Biostrings::chartr("*", "X", translations_combined)

# write fasta of translated sequences
translations_combined %>%
    ape::as.AAbin() %>%
    ape::write.FASTA(., file = "translations.fasta")

# write deduplicated nucleotide sequences
seqs_dss %>% 
    as.character() %>% 
    write_fasta(., "nuc_deduplicated.fasta")
