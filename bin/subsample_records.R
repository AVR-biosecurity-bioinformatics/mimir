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
    "subsample_seeds",
    "subsample_size"
    )
lapply(nf_vars, nf_var_check)

### process variables 

seeds <- readr::read_lines(subsample_seeds) %>% as.integer()
subsample_size <- as.numeric(subsample_size)

# read fasta from file
seqs <- Biostrings::readDNAStringSet(fasta_file)

allowed_ranks <- c("species","genus", "family", "order", "class", "phylum", "kingdom")

### run code

# handle when there are fewer or the same number of sequences as the subsample size
if (length(seqs) <= subsample_size){

    # write all sequences to file
    Biostrings::writeXStringSet(seqs, "subsample0.fasta")

} else {

    # get lineages at every rank for all sequences, and get probability weights based on taxonomic diversity
    lin_total <- 
        names(seqs) %>%
        tibble::as_tibble_col(column_name = "name") %>%
        dplyr::mutate(
            species = stringr::str_extract(name, "(?<=;).+$"),
            genus = stringr::str_extract(species, ".+(?=;)"),
            family = stringr::str_extract(genus, ".+(?=;)"),
            order = stringr::str_extract(family, ".+(?=;)"),
            class = stringr::str_extract(order, ".+(?=;)"),
            phylum = stringr::str_extract(class, ".+(?=;)"),
            kingdom = stringr::str_extract(phylum, ".+(?=;)")
        ) %>%
        # group by each rank's parent lineage, then count distinct taxa per rank
        dplyr::mutate(count_k = dplyr::n_distinct(kingdom)) %>%
        dplyr::mutate(.by = allowed_ranks[7], count_p = dplyr::n_distinct(phylum)) %>%
        dplyr::mutate(.by = allowed_ranks[6], count_c = dplyr::n_distinct(class)) %>%
        dplyr::mutate(.by = allowed_ranks[5], count_o = dplyr::n_distinct(order)) %>%
        dplyr::mutate(.by = allowed_ranks[4], count_f = dplyr::n_distinct(family)) %>%
        dplyr::mutate(.by = allowed_ranks[3], count_g = dplyr::n_distinct(genus)) %>%
        dplyr::mutate(.by = allowed_ranks[2], count_s = dplyr::n_distinct(species)) %>%
        # convert counts to probability per species lineage
        dplyr::mutate(
            weight = (1/count_k) * (1/count_p) * (1/count_c) * (1/count_o) * (1/count_f) * (1/count_g) * (1/count_s),
            # downweight partially classified sequences since they give us less information
            weight = dplyr::if_else(stringr::str_detect(species, ";Unclassified"), weight/10000, weight)
        ) %>%
        dplyr::select(name, weight, species)

    for (i in seeds){
        set.seed(i)
        # randomly group sequences within species into pairs
        lin_ss <- 
            lin_total %>%
            dplyr::slice_sample(prop = 1) %>% # randomise row order
            dplyr::group_by(species) %>%
            dplyr::mutate(pair = paste0(species, ceiling(dplyr::row_number()/2))) %>%
            dplyr::ungroup() %>%
            dplyr::distinct(name, pair, weight)
        # get tibble of just pairs and weights
        lin_pairs <- lin_ss %>% dplyr::distinct(pair, weight)
        # sample as many pairs as required sequences, in the edge case that every selected pair is a singleton
        pairs <- sample(lin_pairs$pair, size = subsample_size, replace = F, prob = lin_pairs$weight)
    
    
        # get sequence names for each pair and only keep the required number
        ss <- lin_ss %>% dplyr::filter(pair %in% pairs) %>% .$name %>% .[1:subsample_size] 
        # subset sequence DSS object
        ss_seqs <- seqs[names(seqs) %in% ss]
        # write sequences to file
        Biostrings::writeXStringSet(ss_seqs, stringr::str_glue("subsample{i}.fasta"))
    }
}



# # get lineages at every rank for all sequences
# lin_total <- 
#     names(seqs) %>%
#     tibble::as_tibble_col(column_name = "name") %>%
#     dplyr::mutate(
#         species = stringr::str_extract(name, "(?<=;).+$"),
#         genus = stringr::str_extract(species, ".+(?=;)"),
#         family = stringr::str_extract(genus, ".+(?=;)"),
#         order = stringr::str_extract(family, ".+(?=;)"),
#         class = stringr::str_extract(order, ".+(?=;)"),
#         phylum = stringr::str_extract(class, ".+(?=;)"),
#         kingdom = stringr::str_extract(phylum, ".+(?=;)")
#     )

# for (i in seeds){
#     set.seed(i)
#     sel_seqs <- character()
#     lin <- lin_total
#     while (length(sel_seqs) < subsample_size){
#         # # set sample size depending on how many sequences remain
#         ### if this is omitted then final subsample size will be subsample_size OR subsample_size + 1 -- but makes code run slightly faster so maybe worth it
#         # if (subsample_size - length(sel_seqs) == 1){
#         #     samp_size <- 1
#         # } else {
#         #     samp_size <- 2
#         # }
#         # sample one kingdom
#         lin_k <- lin[lin$kingdom == sample(lin$kingdom %>% unique(), size = 1),]
#         # sample one phylum
#         lin_p <- lin_k[lin_k$phylum == sample(lin_k$phylum %>% unique(), size = 1),]
#         # sample one class
#         lin_c <- lin_p[lin_p$class == sample(lin_p$class %>% unique(), size = 1),]
#         # sample one order
#         lin_o <- lin_c[lin_c$order == sample(lin_c$order %>% unique(), size = 1),]
#         # sample one family
#         lin_f <- lin_o[lin_o$family == sample(lin_o$family %>% unique(), size = 1),]
#         # sample one genus
#         lin_g <- lin_f[lin_f$genus == sample(lin_f$genus %>% unique(), size = 1),]
#         # sample one species
#         lin_s <- lin_g[lin_g$species == sample(lin_g$species %>% unique(), size = 1),]
#         # sample up to two sequence rows
#         sel <- lin_s %>% dplyr::slice_sample(n = 2) %>% .$name
#         # remove selected sequences from lineage tibble
#         lin <- lin[!lin$name %in% sel,]
#         # update sel_seqs
#         sel_seqs <- c(sel_seqs, sel)
#     }
#     # get final sequence subsample
#     seqs_sub <- seqs[names(seqs) %in% sel_seqs]
#     # write sequences to file
#     Biostrings::writeXStringSet(seqs_sub, stringr::str_glue("subsample{i}.fasta"))
# }