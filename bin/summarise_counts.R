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
    "params_dict",
    "count_genbank",
    "count_bold",
    "count_mito",
    "count_genome",
    "count_internal",
    "count_external",
    "count_input",
    "count_filter_phmm",
    "count_filter_stop",
    "count_remove_exact",
    "count_remove_tax_outliers",
    "count_remove_seq_outliers",
    "count_prune_groups"
    )
lapply(nf_vars, nf_var_check)

### process variables 


### run code

## summarise
counts_summary <- 
    tibble::tibble(
        genbank = count_genbank,
        bold = count_bold,
        mito = count_mito,
        genome = count_genome,
        internal = count_internal,
        external = count_external,
        input = count_input,
        filter_phmm = count_filter_phmm,
        filter_stop = count_filter_stop,
        remove_exact = count_remove_exact,
        remove_tax_outliers = count_remove_tax_outliers,
        remove_seq_outliers = count_remove_seq_outliers,
        prune_groups = count_prune_groups
    ) %>%
    tidyr::pivot_longer(
        cols = genbank:prune_groups,
        names_to = "step",
        values_to = "sequences"
    ) %>%
    dplyr::mutate(
        sequences = as.numeric(sequences),
        step = forcats::fct_relevel(
            step, 
            c(
                "genbank",
                "bold",
                "mito",
                "genome",
                "internal",
                "external",
                "input",
                "filter_phmm",
                "filter_stop",
                "remove_exact",
                "remove_tax_outliers",
                "remove_seq_outliers",
                "prune_groups"
            )
        )
    )

readr::write_csv(counts_summary, "counts_summary.csv")

# plot
ggplot2::ggplot(counts_summary, aes(x = sequences, y = step)) +
  geom_col() +
  geom_text(aes(label = sequences, x = sequences + (max(sequences)* 0.02)), hjust = 0) + # label 5% past end of column 
  scale_y_discrete(limits = rev)

ggsave("counts_summary.pdf", width = 4, height = 6)