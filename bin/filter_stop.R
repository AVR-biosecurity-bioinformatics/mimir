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
# seqs <- readRDS(seqs_file)
seqs <- ape::read.FASTA(fasta_file)

### run code

## run filter only if marker is a coding sequence
if ( coding == "true" ){
    
    ## get genetic code per sequence
    
    # read in ncbi_gencodes from file
    ncbi_gencodes <- readRDS(ncbi_gencodes) 

    allowed_ranks <- c("kingdom","phylum", "class", "order", "family", "genus", "species")

    # get the taxonomic lineage of each sequence as a tibble/df
    lineage <- 
        names(seqs) %>%
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
    if ( length(seqs) != length(genetic_code_v) ){
        stop("Number of sequences in input .fasta does not match number of returned genetic codes")
    }

    # clean up memory
    rm(match_all)
    rm(match_censor_s)
    rm(match_censor_g)
    rm(match_censor_f)
    rm(match_censor_o)
    rm(match_censor_c)
    rm(match_censor_p)
    rm(match_censor_k)
    gc()

    ## do filtering of sequences containing stop codons

    ### NOTE: This assumes all sequences are in the correct orientation ie. 5'-3'.

    # Convert to DNAStringSet
    seqs <- DNAbin2DNAstringset(seqs, remove_gaps=TRUE)

    # get forward reading frames
    frames <- lapply(1:3, function(pos) XVector::subseq(seqs, start=pos))

    ## restructure list to group sequences together (rather than frames)

    # list where each element is a DNAStringSet of the frames of each sequence (at the moment 3 frames hard-coded)
    frames_list <- 
        frames %>% 
        as.data.frame(col.names = c("frame0","frame1","frame2")) %>% 
        tibble::as_tibble(rownames = "name") %>%
        tidyr::pivot_longer(cols = frame0:frame2, names_to = "frame", values_to = "sequence") %>%
        split(., .$name) %>%
        lapply(
            ., 
            function(x) {
                x_df <- as.data.frame(x)
                seqs_DSS <- DNAStringSet(x_df$sequence)
                names(seqs_DSS) <- x$name
                return(seqs_DSS)
            } 
        )

    # check identical length
    if ( length(frames_list) != length(genetic_code_v) ) {
        stop("'frames_list' and 'genetic_code_v' are of different lengths when they should be equal!")
    }

    # translate frames_list using vector of genetic codes
    translated <- 
        mapply( # applies function using combination of two vectors 
            function(y, z) {
                suppressWarnings(
                    Biostrings::translate(
                        x = y, 
                        genetic.code = Biostrings::getGeneticCode(z),
                        if.fuzzy.codon=c("solve", "X")
                    )
                )
            }, 
            frames_list, # y
            genetic_code_v # z
        )

    # check if at least one frame per input sequence contains no stop codons
    stop_codon_retain <- lapply(
        translated, 
        function(a) {
            # get the number of stop codons per translation
            fvec <- c(
                stringr::str_count(as.character(a[1]), "\\*"),
                stringr::str_count(as.character(a[2]), "\\*"),
                stringr::str_count(as.character(a[3]), "\\*")
            )
            # if at least one frame has zero stop codons, retained is true
            if (sum(fvec == 0) >= 1) { 
                retained <- TRUE
            } else {
                retained <- FALSE
            } 
            return(retained)
        }
        ) %>% 
        unlist() %>% 
        unname()

    seqs_filtered <- seqs[stop_codon_retain] # sequences filtered of stop codons


    # ## filter out sequences with stop codons
    # seqs_filtered <- 
    #     codon_filter(
    #         x = seqs, 
    #         genetic_code = params.genetic_code, 
    #         tryrc = TRUE, 
    #         resolve_draws = "majority"
    #     )
} else if ( coding == "false" ) {
    seqs_filtered <- seqs
} else {
    stop("ERROR: 'coding' must be 'true' or 'false'")
}

# write fasta
if ( !is.null(seqs_filtered) && length(seqs_filtered) > 0 ){
    write_fasta(
        seqs_filtered, 
        file = paste0("filter_stop.fasta"), 
        compress = FALSE
        )
} else {
    file.create(paste0("filter_stop.fasta"))
}