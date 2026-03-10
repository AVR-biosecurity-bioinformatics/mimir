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
	"igraph",
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
    "comparators_cfa",
    "cg_csv",
    "seqs_file",
    "counts_file",
	"thresholds_csv",
	"component_group_size"
)
lapply(nf_vars, nf_var_check)

### process variables 

# read in consistent genera csv
cg_all <- readr::read_csv(cg_csv, show_col_types = FALSE)

# read in sequences passing max thresholds
seqs <- Biostrings::readDNAStringSet(seqs_file)

# read in aligned genera .cfa as a list of DSS alignments
comparators_list <- 
	readr::read_lines(comparators_cfa) %>%
	lapply(
		.,
		function(x){
			x_split <- stringr::str_split_1(x, ">>>")
			x_head <- x_split[c(TRUE,FALSE)] %>% stringr::str_remove(., "^>") 
			x_seq <- x_split[c(FALSE,TRUE)]
			names(x_seq) <- x_head
			out <- Biostrings::DNAStringSet(x_seq)
			return(out)
		}
	)

# read in table of redundancy counts (pre-synthetic renaming)
rf_counts <- readr::read_tsv(counts_file, col_names = c("name","n"), show_col_types = FALSE)

# read in thresholds csv
thresholds <- readr::read_csv(thresholds_csv,  show_col_types = FALSE)

# component group size parameters
component_group_size <- as.numeric(component_group_size)

### run code




##### NOTE: select_min_comparators is selecting comparators that are unclassified at the LSR rank, and also only selecting a single comparator per LSR rank

stop()