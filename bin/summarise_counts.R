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

# import counts.csv
counts_raw <- readr::read_csv("counts.csv")

### run code

# summarise
counts_summary <- 
    counts_raw %>%
    dplyr::arrange(order) %>%
    dplyr::mutate(
        process = forcats::fct_reorder(process, order)
    ) %>%
    dplyr::select(process, sequences)

# write .csv
readr::write_csv(counts_summary, "counts_summary.csv")

# plot
counts_plot <- 
    counts_summary %>% 
    ggplot2::ggplot(., aes(x = sequences, y = process)) +
    geom_col() +
    geom_text(
        aes(label = sequences, x = sequences + (max(sequences)* 0.02)), # label 2% past end of column
        hjust = 0, 
        size = 2
    ) +  
    scale_y_discrete(limits = rev) +
    scale_x_continuous(limits = c(0,max(counts_summary$sequences) * 1.2))

# write plot
ggsave("counts_summary.pdf", counts_plot, width = 4, height = 6)