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
    "summary_csv",
    "min_k",
    "max_k"
    )
lapply(nf_vars, nf_var_check)

### process variables 

min_k <- as.numeric(min_k)
max_k <- as.numeric(max_k)

allowed_ranks <- c("species","genus", "family", "order", "class", "phylum", "kingdom")

# read in summary_csv file
summary_tibble <- readr::read_csv(summary_csv)

### run code

# estimate global thresholds from the subsample estimates
summary_summary <- 
    summary_tibble %>%
    dplyr::summarise(
        .by = lsr,
        med = mean(med),
        mad = median(mad)
    ) %>%
    tidyr::pivot_longer(cols = med:mad, names_to = "type", values_to = "value")

# plot the distribution of estimates with the global values
est_plot <- 
    summary_tibble %>%
    tidyr::pivot_longer(cols = med:mad, names_to = "type", values_to = "value") %>%
    ggplot(aes(x = forcats::fct_relevel(lsr, allowed_ranks), y = value)) +
    geom_point(alpha = 0.05, position = position_jitter()) +
    geom_violin(scale = "width", adjust = 2, colour = NA, fill = "lightblue2") +
    geom_point(inherit.aes = F, data = summary_summary,  aes(x = lsr, y = value), colour = "red", shape = 15, size = 2.5) +
    scale_y_continuous(limits = c(NA,NA), name = "value") +
    scale_x_discrete(name = "Lowest shared rank") +
    facet_wrap(.~type, scales = "free")

ggsave("estimation_plot.pdf", est_plot, width = 8, height = 5)

# plot the distribution of the number of pairwise comparisons per LSR across all the subsamples
# Note lower ranks might have much lower sample sizes than higher ranks
n_hist_plot <- 
  summary_tibble %>%
  ggplot(aes(x = n)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(name = "Pairwise comparisons") +
  facet_wrap(.~forcats::fct_relevel(lsr, allowed_ranks), scales = "free", nrow = 2)

ggsave("sample_size_histogram.pdf", n_hist_plot, width = 12, height = 6)

# get thresholds from median, mad and k values
thresholds <- 
    summary_summary %>%
    tidyr::pivot_wider(names_from = "type", values_from = "value") %>%
    dplyr::mutate(
        max = med + max_k*mad, 
        min = med - min_k*mad,
        max = dplyr::case_when(
            is.na(max) ~ 0,
            max > 1 ~ 1,
            max < 0 ~ 0,
            .default = max
        ),
        min = dplyr::case_when(
            is.na(min) ~ 0,
            min > 1 ~ 1,
            min < 0 ~ 0,
            .default = min
        )
    )

readr::write_csv(thresholds, "thresholds.csv")