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
    "ggplot2",
    "httr",
    "kmer",
    "magrittr",
    "methods",
    "openssl",
    "parallel",
    "patchwork",
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
    "params_dict"
    )
lapply(nf_vars, nf_var_check)

### process variables 
sources_fates <- readr::read_csv("sources_fates.csv")

# columns of sources
cols_sources <- c(
    genbank = 0,
    bold = 0,
    mito = 0, 
    genome = 0, 
    internal = 0
)

# order of processes
order_sources <- c(
    "genbank",
    "bold",
    "mito",
    "genome",
    "internal"
)

# columns of fates
cols_fates <- c(
    filter_unclassified = 0,
    filter_phmm_full = 0,
    filter_phmm_trimmed = 0, 
    filter_duplicates = 0, 
    filter_ambiguous = 0, 
    filter_tax_outliers = 0, 
    filter_seq_outliers = 0, 
    filter_redundant = 0, 
    final_database = 0
)

# order of fates
order_fates <- c(
    "filter_unclassified",
    "filter_phmm_full",
    "filter_phmm_trimmed",
    "filter_duplicates",
    "filter_ambiguous",
    "filter_tax_outliers",
    "filter_seq_outliers",
    "filter_redundant",
    "final_database"
)

### run code 

# extract metadata from each sequence name 
sf_meta <- 
    sources_fates %>%
    tidyr::separate(name, into = c("id", "kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", remove = FALSE) %>%
    tidyr::separate(id, into = c("accession", "taxid"), sep = "\\|", remove = TRUE) %>%
    dplyr::mutate(
        source = forcats::fct_relevel(source, order_sources), 
        fate = forcats::fct_relevel(fate, order_fates)
    ) %>%
    dplyr::arrange(kingdom, phylum, class, order, family, genus, species, source, fate) 

# save metadata
readr::write_csv(sf_meta, "sf_meta.csv")


# get total number of sequences
n_seqs_input <- nrow(sf_meta)
n_seqs_final <- sf_meta %>% dplyr::filter(fate == "final_database") %>% nrow()

# summarise number of sequences from each source
sf_meta %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = source, values_from = n, values_fill = 0) %>%
    tibble::add_column(!!!cols_sources[!names(cols_sources) %in% names(.)]) %>%
    tidyr::pivot_longer(cols = names(cols_sources), names_to = "source", values_to = "sequences")

# summarise number of sequences with each fate
# using accumulate
fate_summary <- 
    sf_meta %>%
    dplyr::group_by(fate) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>% 
    tidyr::pivot_wider(names_from = fate, values_from = n, values_fill = 0) %>%
    tibble::add_column(!!!cols_fates[!names(cols_fates) %in% names(.)]) %>%
    tidyr::pivot_longer(cols = names(cols_fates), names_to = "fate", values_to = "removed") %>%
    dplyr::mutate(
        input = purrr::accumulate(
        .x = removed, 
        .init = n_seqs_input,
        .f = ~ .x - .y
        )[1:n()],
        retained = input - removed
    ) %>%
    dplyr::mutate(
        fate = forcats::fct_relevel(fate, order_fates)
    ) %>%
    tibble::add_row(
        fate = "input",
        removed = 0,
        retained = n_seqs_input,
        .before = 1
    ) %>%
    dplyr::mutate(fate = forcats::fct_relevel(fate, rev(c("input", order_fates)))) %>%
    dplyr::mutate(
        removed = dplyr::case_when(fate == "final_database" ~ 0, .default = removed),
        retained = dplyr::case_when(fate == "final_database" ~ n_seqs_final, .default = retained)
    )

## plot fates
# plot retained
f_ret_plot <- 
    ggplot(fate_summary, aes(y = fate, x = retained)) +
    geom_col(fill = "skyblue1") +
    geom_text(aes(label = retained, x = retained + (0.05*n_seqs_input)), hjust = 0, size = 2 ) +
    scale_x_continuous(limits = c(0,n_seqs_input*1.2)) +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

# plot removed
f_rem_plot <- 
    ggplot(fate_summary, aes(y = fate, x = removed)) +
    geom_col(fill = "lightcoral") +
    geom_text(aes(label = dplyr::if_else(fate %in% c("input", "final_database"), NA, removed), x = removed + (0.05*n_seqs_input)), hjust = 1, size = 2 ) +
    scale_x_reverse(limits = c(n_seqs_input*1.2,0)) +
    scale_y_discrete(name = "Step")

# plot combined
fates_plot <- (f_rem_plot + f_ret_plot)

ggsave("fates_plot.pdf", fates_plot, height = 4, width = 6)


## taxonomic summary
taxa_summary <- 
    sf_meta %>%
    dplyr::group_by(fate) %>%
    dplyr::summarise_all(n_distinct) %>%
    dplyr::select(
        fate,
        kingdom, 
        phylum,
        class,
        order,
        family,
        genus,
        species
    ) %>%
    dplyr::arrange(fate)

readr::write_csv(taxa_summary, "taxa_summary.csv")