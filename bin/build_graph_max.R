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
    "flagged_genera_file",
    "seqs_file",
    "counts_file",
	"component_group_size"
)
lapply(nf_vars, nf_var_check)

### process variables 

flagged_genera <- readr::read_csv(flagged_genera_file, show_col_types = FALSE)

seqs <- Biostrings::readDNAStringSet(seqs_file)

counts <- readr::read_tsv(counts_file, col_names = c("name","n"), show_col_types = FALSE)

component_group_size <- as.numeric(component_group_size)

### run code

# make tibble of genus lineage of each sequence record
seqs_genus <- 
	seqs %>%
	names() %>%
	tibble::as_tibble_col(column_name = "name") %>%
    tidyr::separate(col = name, into = c("seqid", "lineage_string"), sep = ";", extra = "merge", remove = F) %>%
    dplyr::mutate(
        genus = stringr::str_extract(lineage_string, "^.+(?=;)")
    ) %>%
	dplyr::select(name, genus)

## make sure counts sequence names match new names for synthetic genera
# LCR family or above
counts_lcrf <- 
	counts %>%
	dplyr::filter(stringr::str_detect(name, ";Unclassified;[^;]+$"))

# LCR genus or below
count_lcrg <- 
	counts %>%
	dplyr::filter(stringr::str_detect(name, ";Unclassified;[^;]+$", negate = T))

# make new counts tibble
counts_new <- 
	names(seqs)[stringr::str_detect(names(seqs), ";Unclassified\\d+;[^;]+$")] %>%
	tibble::as_tibble_col(column_name = "new") %>%
	dplyr::mutate(
		name = stringr::str_replace(new, ";Unclassified\\d+;(?=[^;]+$)", ";Unclassified;")
	) %>%
	dplyr::left_join(., counts_lcrf, by = "name") %>%
	dplyr::select(name = new, n) %>%
	dplyr::bind_rows(., count_lcrg)

# summarise counts by genus
genus_summary <- 
	counts_new %>%
	dplyr::mutate(genus = stringr::str_remove(name, "^[^;]+;") %>% stringr::str_remove(., ";[^;]+$")) %>%
	dplyr::summarise(
		.by = genus, 
		rep_n = sum(n),
		act_n = n()
	) %>%
	dplyr::arrange(desc(act_n))

# build graph of max violations and export components
max_graph <- 
	flagged_genera %>%
	dplyr::distinct() %>%
	igraph::graph_from_data_frame(
        ., 
        directed = F
    ) %>%
    igraph::simplify() 

# determine each component of graph
max_components <- 
    max_graph %>%
    igraph::components() %>%
    .$membership %>%
    tibble::enframe(name = "genus", value = "component") %>%
    dplyr::arrange(component, genus) %>%
	dplyr::left_join(., genus_summary, by = "genus")

# group small components together so groups with more than 1 member have max 1000 sequences
# this is so we don't align more than 1000 sequences together unless we have to
component_groups <- 
	max_components %>%
	dplyr::summarise(
		.by = component, 
		rep_n = sum(rep_n),
		act_n = sum(act_n),
		genus_n = n()
	) %>%
	dplyr::arrange(act_n) %>%
	dplyr::mutate(
		# idea from here: https://stackoverflow.com/questions/34531568/conditional-cumsum-with-reset
		cs = base::Reduce(\(x, y) if (x + y <= 1000) x + y else y, x = act_n, accumulate = TRUE),
		group = cumsum(act_n == cs)
	) %>%
	dplyr::select(component, group) 

mc_split <- 
	max_components %>%
	dplyr::left_join(., component_groups, by = "component") %>%
	dplyr::arrange(group, component) %>%
	split(., .$group) 

# for each group of components of genera, get the associated sequence records and write to file
lapply(
	seq_along(mc_split),
	function(y){
		x <- mc_split[[y]]
		x_seq_names <- 
			seqs_genus %>%
			dplyr::filter(genus %in% x$genus) %>%
			dplyr::pull(name) 
		x_seqs <- seqs[names(seqs) %in% x_seq_names]
		Biostrings::writeXStringSet(x_seqs, stringr::str_glue("component_group{y}.fasta"), width = 9999)
		message(stringr::str_glue("{length(x_seqs)} sequences in component group {y}"))
	}
)