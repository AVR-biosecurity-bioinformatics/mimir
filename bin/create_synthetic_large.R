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
    "fasta_files",
    "clusters_tsv",
    "rf_counts_tsv",
    "small_reps_txt"
)
lapply(nf_vars, nf_var_check)

### process variables 

# read in list of sequences
seqs_list <-
    stringr::str_extract_all(fasta_files, pattern = "[^\\s,\\[\\]]+") %>% 
    unlist() %>%
    lapply(., Biostrings::readDNAStringSet) 

# read in clusters tsv
clusters <- readr::read_tsv(clusters_tsv, col_names = c("rep","name"), show_col_types = FALSE)

# read in rf_counts csv
rf_counts <- readr::read_tsv(rf_counts_tsv, col_names = c("name","n"), show_col_types = FALSE)

# read in rep names from small synthetic genera
small_reps <- readr::read_lines(small_reps_txt)

### run code

# get the max cluster ID number (from small clusters) for each genus 
genera_renamed <- 
	small_reps %>%
	tibble::as_tibble_col(column_name = "name") %>%
	dplyr::mutate(
		genus_new = stringr::str_extract(name, "(?<=;).+$") %>% stringr::str_extract(., ".+(?=;)"),
		genus_old = stringr::str_remove(genus_new, "\\d+$"),
		cluster_number = stringr::str_extract(genus_new, "\\d+$") %>% as.numeric()
	) %>%
	dplyr::mutate(
		.by = genus_old,
		max = max(cluster_number) %>% as.integer()
	) %>%
	dplyr::distinct(genus_old, max)

# create synthetic genera for large clusters
clusters_large_renamed <- 
	clusters %>%
	dplyr::mutate(
		.by = rep,
		# get old genus name
		genus_old = stringr::str_extract(name, "(?<=;).+$") %>% stringr::str_extract(., ".+(?=;)")
	) %>%
	dplyr::arrange(genus_old) %>%
	# add max ID for each genus
	dplyr::left_join(., genera_renamed, by = "genus_old") %>%
	dplyr::mutate(
		.by = rep,
		# use max to make sure new cluster ID is always higher than those in small clusters
		genus_new = paste0(genus_old, dplyr::cur_group_id() + max),
		# create new sequence name using new genus
		name_new = stringr::str_replace(name, paste0(";",genus_old), paste0(";",genus_new))
	)

# use new names to rename representative counts and export just changed lines
rf_counts_renamed <- 
	clusters_large_renamed %>%
	dplyr::left_join(., rf_counts, by = "name") %>%
	dplyr::select(name = name_new, n)

readr::write_tsv(rf_counts_renamed, "counts_renamed.tsv", col_names = FALSE)

# rename large clusters as synthetic genera and create DSS object for export
seqs_synthetic <- 
	seqs_list %>% 
	# collect all sequences into a single file
	do.call(c, .) %>%
	.[names(.) %in% clusters_large_renamed$name] %>% 
	as.character() %>% 
	tibble::enframe() %>%
	dplyr::left_join(., clusters_large_renamed, by = "name") %>%
	dplyr::select(name_new, value) %>%
	tibble::deframe() %>%
	Biostrings::DNAStringSet()

if(length(seqs_synthetic) > 0){
	Biostrings::writeXStringSet(seqs_synthetic, "synthetic_genera.fasta")
} else {
	file.create("synthetic_genera.fasta")
}

# output representative sequences names for each cluster
rep_names <- 
	clusters_large_renamed %>%
	dplyr::distinct(rep) %>%
	dplyr::select(name = rep) %>%
	dplyr::left_join(., clusters_large_renamed, by = "name") %>%
	dplyr::pull(name_new)

readr::write_lines(rep_names, "reps.txt")