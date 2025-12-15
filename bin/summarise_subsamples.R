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
    "fasta_files"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# read in list of sequences
fasta_list <- # convert Groovy to R list format
    stringr::str_extract_all(fasta_files, pattern = "[^\\s,\\[\\]]+") %>% unlist()

seqs_list <- lapply(fasta_list, Biostrings::readDNAStringSet) 

allowed_ranks <- c("kingdom","phylum", "class", "order", "family", "genus", "species")
root_ranks <- c("root", allowed_ranks)

joint_summary <- 
    lapply(
        1:length(seqs_list),
        function(i){
            x <- seqs_list[[i]]
            x_d <- DECIPHER::DistanceMatrix(x, verbose = F)
            # remove redundant pairwise comparisons, including self-comparisons
            x_d[lower.tri(x_d, diag = TRUE)] <- NA
            # get pairwise comparisons from matrix
            x_p <- 
                lapply(
                    1:nrow(x_d),
                    function(y){
                        # get taxonomy of focal sequence
                        idy <- rownames(x_d[y, , drop = F])
                        id_tax <- 
                            stringr::str_split(idy, ";", n = 8, simplify = F) %>% 
                            unlist %>% 
                            stringr::str_replace_all(., "^Unclassified$", "UNCLASSIFIED") %>% 
                            .[2:8] %>%
                            c("Root", .)
                        id_r <- "Root"
                        id_k <- base::paste(id_tax[1:2], collapse = ";")
                        id_p <- base::paste(id_tax[1:3], collapse = ";")
                        id_c <- base::paste(id_tax[1:4], collapse = ";")
                        id_o <- base::paste(id_tax[1:5], collapse = ";")
                        id_f <- base::paste(id_tax[1:6], collapse = ";")
                        id_g <- base::paste(id_tax[1:7], collapse = ";")
                        id_s <- base::paste(id_tax[1:8], collapse = ";")
                        id_classified <- stringr::str_detect(id_tax, "UNCLASSIFIED", negate = T)
                        # get lowest classified rank for focal sequence (ie. max index)
                        id_lcr <- root_ranks[max(match(root_ranks[id_classified], root_ranks))] 
                        id_lcr_i <- which(root_ranks == id_lcr)
                        
                        init <-
                            x_d[y, ] %>% 
                            tibble::enframe() %>%
                            tidyr::drop_na() 
                        
                        if (nrow(init) == 0){
                            y_parsed <- 
                                tibble::tibble(
                                    fid = numeric(),
                                    lsr = character(),
                                    lcr = character()
                                )
                        } else {
                            y_parsed <- 
                                init %>%
                                dplyr::mutate(
                                    fid = 1 - value, # convert distance to identity
                                    species = paste("Root", stringr::str_extract(name, "(?<=;).+$"), sep = ";"),
                                    genus = stringr::str_extract(species, ".+(?=;)"),
                                    family = stringr::str_extract(genus, ".+(?=;)"),
                                    order = stringr::str_extract(family, ".+(?=;)"),
                                    class = stringr::str_extract(order, ".+(?=;)"),
                                    phylum = stringr::str_extract(class, ".+(?=;)"),
                                    kingdom = stringr::str_extract(phylum, ".+(?=;)"),
                                    # lowest shared rank (discounting unclassified ranks in either sequence)
                                    lsr = dplyr::case_when(
                                        id_s == species ~ "species",
                                        id_g == genus ~ "genus",
                                        id_f == family ~ "family",
                                        id_o == order ~ "order",
                                        id_c == class ~ "class",
                                        id_p == phylum ~ "phylum",
                                        id_k == kingdom ~ "kingdom",
                                        .default = NA
                                    ),
                                    # lowest classified rank for the target sequence
                                    lcr = dplyr::case_when(
                                        stringr::str_detect(species, "Unclassified", negate = T) ~ "species",
                                        stringr::str_detect(genus, "Unclassified", negate = T) ~ "genus",
                                        stringr::str_detect(family, "Unclassified", negate = T) ~ "family",
                                        stringr::str_detect(order, "Unclassified", negate = T) ~ "order",
                                        stringr::str_detect(class, "Unclassified", negate = T) ~ "class",
                                        stringr::str_detect(phylum, "Unclassified", negate = T) ~ "phylum",
                                        stringr::str_detect(kingdom, "Unclassified", negate = T) ~ "kingdom",
                                        .default = "root"
                                    )
                                ) %>%
                                dplyr::rowwise() %>%
                                dplyr::mutate(
                                    # get the lowest rank that is classified for both sequences
                                    lcr = root_ranks[min(id_lcr_i, which(root_ranks == lcr))]
                                ) %>%
                                dplyr::ungroup() %>%
                                dplyr::select(fid, lsr, lcr)
                        }
                        return(y_parsed)
                    }
                ) %>%
                dplyr::bind_rows()
                
                # species LSR
                lsr_s <- 
                    x_p %>%
                    dplyr::filter(
                        lsr == "species" & lcr %in% c("species")
                    ) %>%  
                    dplyr::summarise(
                        lsr = "species",
                        n = n(),
                        med = median(fid),
                        mad = mad(fid)
                    )
                # genus LSR
                lsr_g <- 
                    x_p %>%
                    dplyr::filter(
                        lsr == "genus" & lcr %in% c("species")
                    ) %>%
                    dplyr::summarise(
                        lsr = "genus",
                        n = n(),
                        med = median(fid),
                        mad = mad(fid)
                    )
                # family LSR
                lsr_f <- 
                    x_p %>%
                    dplyr::filter(
                        lsr == "family" & lcr %in% c("species","genus")
                    ) %>%
                    dplyr::summarise(
                        lsr = "family",
                        n = n(),
                        med = median(fid),
                        mad = mad(fid)
                    )
                # order LSR
                lsr_o <- 
                    x_p %>%
                    dplyr::filter(
                        lsr == "order" & lcr %in% c("species","genus", "family")
                    ) %>%
                    dplyr::summarise(
                        lsr = "order",
                        n = n(),
                        med = median(fid),
                        mad = mad(fid)
                    )
                # class LSR
                lsr_c <- 
                    x_p %>%
                    dplyr::filter(
                        lsr == "class" & lcr %in% c("species","genus", "family", "order")
                    ) %>%
                    dplyr::summarise(
                        lsr = "class",
                        n = n(),
                        med = median(fid),
                        mad = mad(fid)
                    )
                # phylum LSR
                lsr_p <- 
                    x_p %>%
                    dplyr::filter(
                        lsr == "phylum" & lcr %in% c("species","genus", "family", "order", "class")
                    ) %>%
                    dplyr::summarise(
                        lsr = "phylum",
                        n = n(),
                        med = median(fid),
                        mad = mad(fid)
                    )
                # kingdom LSR
                lsr_k <- 
                    x_p %>%
                    dplyr::filter(
                        lsr == "kingdom" & lcr %in% c("species","genus", "family", "order", "class", "phylum")
                    ) %>%
                    dplyr::summarise(
                        lsr = "kingdom",
                        n = n(),
                        med = median(fid),
                        mad = mad(fid)
                    )
                # combine together
                lsr_all <- 
                    dplyr::bind_rows(
                        lsr_s, 
                        lsr_g, 
                        lsr_f, 
                        lsr_o,
                        lsr_c,
                        lsr_p,
                        lsr_k
                    )

                gc()

                message(stringr::str_glue("Finished file {i} of {length(seqs_list)}"))

                return(lsr_all)
        }
    ) %>%
    dplyr::bind_rows()

# write summary of results to file
readr::write_csv(joint_summary, "joint_summary.csv")