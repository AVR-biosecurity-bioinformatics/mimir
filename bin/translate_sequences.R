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
seqs <- Biostrings::readDNAStringSet(fasta_file)

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

# convert "Unclassified" ranks to "UNCLASSIFIED" in gencodes tibble to prevent exact matching of lower ranks
ncbi_gencodes <- 
    ncbi_gencodes %>%
    dplyr::mutate(
        dplyr::across(
            tidyselect::all_of(allowed_ranks), 
            ~dplyr::if_else(. == "Unclassified", "UNCLASSIFIED", .)
        )
    )

### run code

# deduplicate sequences
seqs_dss <- 
    seqs %>% 
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
    tidyr::separate(col=id, into=c("seqid", "taxid"), sep="\\|", extra="merge") %>%
    dplyr::mutate(seq_order = dplyr::row_number())


# split sequences into those that have an NCBI taxid and those that don't
lineage.taxid_ncbi <- 
    lineage %>%
    dplyr::filter(stringr::str_detect(taxid, "^NCBI:", negate = F))

lineage.taxid_other <- 
    lineage %>%
    dplyr::filter(stringr::str_detect(taxid, "^NCBI:", negate = T))

## for those with NCBI taxid, match exactly using taxid
gc.taxid_ncbi.pre <- 
    lineage.taxid_ncbi %>%
    dplyr::mutate(tax_id = stringr::str_remove(taxid, "^NCBI:") %>% as.integer) %>%
    dplyr::left_join(
        ., 
        ncbi_gencodes %>% 
            dplyr::select(-c(kingdom, phylum, class, order, family, genus, species)), 
        by = "tax_id"
    ) 

# take those that didn't match and combine with non-NCBI
lineage.taxid_other.comb <- 
    gc.taxid_ncbi.pre %>%
    dplyr::filter(is.na(gencode)) %>%
    dplyr::select(seqid, taxid, kingdom, phylum, class, order, family, genus, species, seq_order) %>%
    dplyr::bind_rows(., lineage.taxid_other)

# format matched ncbi
gc.taxid_ncbi <- 
    gc.taxid_ncbi.pre %>%
    dplyr::filter(!is.na(gencode)) %>%
    dplyr::select(seqid, taxid, seq_order, gencode, mtgencode, plgencode, hdgencode)

## for those without NCBI taxid (or didn't match gencodes tibble), match in a different way
# rename sequence ranks to distinguish them from ncbi_gencode columns
lineage.taxid_other.rn <- 
    lineage.taxid_other.comb %>%
    dplyr::rename_with(
        .fn = ~ stringr::str_replace(.x, "$", "_seq"),
        .cols = c(kingdom, phylum, class, order, family, genus, species)
    )

# join sequences to gencodes by each of species, genus and family (individually), keeping the best match across all ranks
species_bm <- 
    ncbi_gencodes %>%
    dplyr::filter(rank == "species") %>%
    dplyr::left_join(., lineage.taxid_other.rn, by = dplyr::join_by("species" == "species_seq"), keep = TRUE) %>% 
    dplyr::filter(!is.na(seqid)) %>%
    # calculate number of rank matches
    dplyr::mutate(
        rank_matches = 
        (kingdom == kingdom_seq) + 
        (phylum == phylum_seq) +
        (class == class_seq) +
        (order == order_seq) +
        (family == family_seq) +
        (genus == genus_seq) +
        (species == species_seq)
    ) %>% 
    # keep best match per sequence
    dplyr::group_by(seqid) %>% 
    dplyr::arrange(desc(rank_matches)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

message("Created species_bm")
gc()

genus_bm <- 
    ncbi_gencodes %>%
    dplyr::group_by(kingdom, phylum, class, order, family, genus, gencode, mtgencode, plgencode, hdgencode) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(., lineage.taxid_other.rn, by = dplyr::join_by("genus" == "genus_seq"), keep = TRUE) %>% 
    dplyr::filter(!is.na(seqid)) %>%
    # calculate number of rank matches
    dplyr::mutate(
        rank_matches = 
        (kingdom == kingdom_seq) + 
        (phylum == phylum_seq) +
        (class == class_seq) +
        (order == order_seq) +
        (family == family_seq) +
        (genus == genus_seq) +
        (species == species_seq)
    ) %>% 
    # keep best match per sequence
    dplyr::group_by(seqid) %>% 
    dplyr::arrange(desc(rank_matches)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() 

message("Created genus_bm")
gc()

family_bm <- 
    ncbi_gencodes %>%
    dplyr::group_by(kingdom, phylum, class, order, family, gencode, mtgencode, plgencode, hdgencode) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(., lineage.taxid_other.rn, by = dplyr::join_by("family" == "family_seq"), keep = TRUE) %>% 
    dplyr::filter(!is.na(seqid)) %>%
    # calculate number of rank matches
    dplyr::mutate(
        rank_matches = 
        (kingdom == kingdom_seq) + 
        (phylum == phylum_seq) +
        (class == class_seq) +
        (order == order_seq) +
        (family == family_seq) +
        (genus == genus_seq) +
        (species == species_seq)
    ) %>% 
    # keep best match per sequence
    dplyr::group_by(seqid) %>% 
    dplyr::arrange(desc(rank_matches)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() 

message("Created family_bm")
gc()

# join sequences to non-redundant gencodes by each of order, class, phylum and kingdom, keeping the best match across higher ranks
# (removing redundancy is done to stop combinatorial explosion destroying memory)
order_bm <- 
    ncbi_gencodes %>%
    dplyr::group_by(kingdom, phylum, class, order, gencode, mtgencode, plgencode, hdgencode) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(., lineage.taxid_other.rn, by = dplyr::join_by("order" == "order_seq"), keep = TRUE) %>% 
    dplyr::filter(!is.na(seqid)) %>%
    # calculate number of rank matches
    dplyr::mutate(
        rank_matches = 
        (kingdom == kingdom_seq) + 
        (phylum == phylum_seq) +
        (class == class_seq) +
        (order == order_seq) +
        (family == family_seq) +
        (genus == genus_seq) +
        (species == species_seq)
    ) %>% 
    dplyr::group_by(seqid) %>% 
    dplyr::arrange(desc(rank_matches)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() 

message("Created order_bm")
gc()

class_bm <- 
    ncbi_gencodes %>%
    dplyr::group_by(kingdom, phylum, class, gencode, mtgencode, plgencode, hdgencode) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(., lineage.taxid_other.rn, by = dplyr::join_by("class" == "class_seq"), keep = TRUE) %>% 
    dplyr::filter(!is.na(seqid)) %>%
    # calculate number of rank matches
    dplyr::mutate(
        rank_matches = 
        (kingdom == kingdom_seq) + 
        (phylum == phylum_seq) +
        (class == class_seq) +
        (order == order_seq) +
        (family == family_seq) +
        (genus == genus_seq) +
        (species == species_seq)
    ) %>% 
    dplyr::group_by(seqid) %>% 
    dplyr::arrange(desc(rank_matches)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() 

message("Created class_bm")
gc()

phylum_bm <- 
    ncbi_gencodes %>%
    dplyr::group_by(kingdom, phylum, gencode, mtgencode, plgencode, hdgencode) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(., lineage.taxid_other.rn, by = dplyr::join_by("phylum" == "phylum_seq"), keep = TRUE) %>% 
    dplyr::filter(!is.na(seqid)) %>%
    # calculate number of rank matches
    dplyr::mutate(
        rank_matches = 
        (kingdom == kingdom_seq) + 
        (phylum == phylum_seq) +
        (class == class_seq) +
        (order == order_seq) +
        (family == family_seq) +
        (genus == genus_seq) +
        (species == species_seq)
    ) %>% 
    dplyr::group_by(seqid) %>% 
    dplyr::arrange(desc(rank_matches)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() 

message("Created phylum_bm")
gc()

kingdom_bm <- 
    ncbi_gencodes %>%
    dplyr::group_by(kingdom, gencode, mtgencode, plgencode, hdgencode) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(., lineage.taxid_other.rn, by = dplyr::join_by("kingdom" == "kingdom_seq"), keep = TRUE) %>% 
    dplyr::filter(!is.na(seqid)) %>%
    # calculate number of rank matches
    dplyr::mutate(
        rank_matches = 
        (kingdom == kingdom_seq) + 
        (phylum == phylum_seq) +
        (class == class_seq) +
        (order == order_seq) +
        (family == family_seq) +
        (genus == genus_seq) +
        (species == species_seq)
    ) %>% 
    dplyr::group_by(seqid) %>% 
    dplyr::arrange(desc(rank_matches)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() 

message("Created kingdom_bm")
gc()

# combine all the *_bm tibbles
lineage.taxid_other.bm <- 
    dplyr::bind_rows(
        species_bm,
        genus_bm,
        family_bm,
        order_bm,
        class_bm,
        phylum_bm,
        kingdom_bm
    ) %>%
    dplyr::group_by(seqid) %>%
    dplyr::arrange(desc(rank_matches)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

rm(species_bm, genus_bm, family_bm, order_bm, class_bm, phylum_bm, kingdom_bm)

# get sequence IDs with rank_matches of <2
bad_matches <- 
    lineage.taxid_other.bm %>%
    dplyr::filter(rank_matches < 2)

if ( nrow(bad_matches) > 0 ){
    message(paste0("Some sequences have unknown genetic codes (max 10 printed here):"))
    message(cat("",head(bad_matches$seqid,10), sep = "\n"))
    stop("Stopping pipeline.")
}

# extract matches from non-NCBI matching
gc.taxid_other <- 
    lineage.taxid_other.bm %>%
    dplyr::select(seqid, taxid, seq_order, gencode, mtgencode, plgencode, hdgencode)

# check all seqids are present
if ( !setequal(gc.taxid_other$seqid, lineage.taxid_other.comb$seqid) ){
    stop(paste0("Some seqids that didn't have NCBI gencodes are missing from the genetic code matching tibble!"))
}

# combine into one tibble in original sequence order
gc.combined <- 
    dplyr::bind_rows(gc.taxid_ncbi, gc.taxid_other) %>%
    dplyr::arrange(seq_order)

# check all seqids are present
if ( !setequal(gc.combined$seqid, lineage$seqid) ){
    paste0("Some seqids are missing from the combined genetic code matching tibble!")
}

# get vector of genetic codes for each sequence depending on marker type
genetic_code_v <- 
    switch(
        marker_type,
        nuclear = gc.combined$gencode,
        mitochondrial = gc.combined$mtgencode,
        plastid = gc.combined$plgencode,
        hydrogen = gc.combined$hdgencode
    )

# check all values of vector are not NA
if (any(is.na(genetic_code_v))){
    stop(paste0("Some genetic codes are undefined!"))
}

# check length of seqs is the same as genetic_code_v
if ( length(seqs_dss) != length(genetic_code_v) ){
    stop("Number of sequences in input .fasta does not match number of returned genetic codes!")
}




# clean up memory
rm(lineage)
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
