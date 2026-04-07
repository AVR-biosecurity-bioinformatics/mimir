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
    "synthetic_max_size"
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

synthetic_max_size <- as.numeric(synthetic_max_size)

### run code

# determine clusters and their sizes
clusters_size <- 
	clusters %>%
	dplyr::mutate(
		.by = rep,
		size = n()
	) 

clusters_small <- clusters_size %>% dplyr::filter(size <= synthetic_max_size)
clusters_large <- clusters_size %>% dplyr::filter(size > synthetic_max_size)

## small clusters

# for small clusters, rename members and output combined .fasta
clusters_small_renamed <- 
	clusters_small %>%
	dplyr::mutate(
		.by = rep,
		# get old genus name
		genus_old = stringr::str_extract(name, "(?<=;).+$") %>% stringr::str_extract(., ".+(?=;)"),
		# make new genus name
		genus_new = paste0(genus_old, dplyr::cur_group_id()),
		# create new sequence name using new genus
		name_new = stringr::str_replace(name, paste0(";",genus_old), paste0(";",genus_new))
	)

# use new names to rename representative counts and export just changed lines
rf_counts_renamed <- 
	clusters_small_renamed %>%
	dplyr::left_join(., rf_counts, by = "name") %>%
	dplyr::select(name = name_new, n)

readr::write_tsv(rf_counts_renamed, "counts_renamed.tsv", col_names = FALSE)

# rename small clusters as synthetic genera and create DSS object for export
seqs_synthetic <- 
	seqs_list %>% 
	# collect all sequences into a single file
	do.call(c, .) %>%
	.[names(.) %in% clusters_small$name] %>% 
	as.character() %>% 
	tibble::enframe() %>%
	dplyr::left_join(., clusters_small_renamed, by = "name") %>%
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
	clusters_small %>%
	dplyr::distinct(rep) %>%
	dplyr::select(name = rep) %>%
	dplyr::left_join(., clusters_small_renamed, by = "name") %>%
	dplyr::pull(name_new)

readr::write_lines(rep_names, "reps.txt")

## large clusters

if(length(clusters_large$name) > 0){
	# for each input .fasta, get the sequences that were found in large clusters and then write to file
	purrr::imap(
		seqs_list,
		function(x, idx){
			seqs_out <- x[names(x) %in% clusters_large$name] 
			if (length(seqs_out) > 0){
				Biostrings::writeXStringSet(
					seqs_out,
					paste0("large",idx,".fasta")
				)
				message(stringr::str_glue("Wrote {length(seqs_out)} sequences to 'large{idx}.fasta'"))
			}
		}
	)
} else {
	# create empty dummy output
	file.create("large0.fasta")
}


