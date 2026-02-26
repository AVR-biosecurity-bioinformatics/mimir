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
    "intragenus_csv",
    "seqs_all_fasta",
    "seqs_max_fasta",
    "aligned_genera_cfa",
    "synthetic_reps_txt",
    "rf_counts_tsv",
    "con_min_n",
    "con_min_prop"
)
lapply(nf_vars, nf_var_check)

### process variables 

# read in intragenus csv
intragenus <- readr::read_csv(intragenus_csv, show_col_types = FALSE)

# read in sequences that passed intragenus filtering as names
seqs_all_names <- Biostrings::readDNAStringSet(seqs_all_fasta) %>% names()

# read in sequences passing max thresholds
seqs_max <- Biostrings::readDNAStringSet(seqs_max_fasta)

# read in aligned genera .cfa as a list of DSS alignments
aligned_genera_list <- 
	readr::read_lines(aligned_genera_cfa) %>%
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

# read in vector of synthetic genus representatives
synthetic_reps <- readr::read_lines(synthetic_reps_txt)

# read in table of redundancy counts (pre-synthetic renaming)
rf_counts <- readr::read_tsv(rf_counts_tsv, col_names = c("name","n"), show_col_types = FALSE)

# consensus parameters
con_min_n <- as.numeric(con_min_n)
con_min_prop <- as.numeric(con_min_prop)

### run code

allowed_ranks <- c("kingdom","phylum", "class", "order", "family", "genus", "species")
root_ranks <- c("root", allowed_ranks)

set.seed(1)

## make sure counts sequence names match new names for synthetic genera
# LCR family or above
counts_lcrf <- 
	rf_counts %>%
	dplyr::filter(stringr::str_detect(name, ";Unclassified;[^;]+$"))

# LCR genus or below
count_lcrg <- 
	rf_counts %>%
	dplyr::filter(stringr::str_detect(name, ";Unclassified;[^;]+$", negate = T))

# make new counts tibble using seqs_all_names
counts_new <- 
	seqs_all_names[stringr::str_detect(seqs_all_names, ";Unclassified\\d+;[^;]+$")] %>%
	tibble::as_tibble_col(column_name = "new") %>%
	dplyr::mutate(
		name = stringr::str_replace(new, ";Unclassified\\d+;(?=[^;]+$)", ";Unclassified;")
	) %>%
	dplyr::left_join(., counts_lcrf, by = "name") %>%
	dplyr::select(name = new, n) %>%
	dplyr::bind_rows(., count_lcrg)


### redetermine intragenus consistency by removing sequences not retained after max outlier detection

# names of sequences removed at max outlier detection
max_removed <- seqs_all_names[!seqs_all_names %in% names(seqs_max)]

## remove max outliers from intragenus results and then recalculate split genera
intragenus_new <- 
	intragenus %>%
	dplyr::filter(threshold == "genus_min") %>%
	dplyr::filter(!name %in% max_removed) %>% 
	# keep only needed columns
	dplyr::select(name, n, taxon, cluster, threshold) %>%
	# recalculate results
	dplyr::mutate(
		.by = c(taxon),
		# is taxon split?
		split = dplyr::if_else(length(unique(cluster)) == 1, FALSE, TRUE),
		# records in genus
		taxon_n = sum(n)
	) %>% 
	dplyr::mutate(
		.by = c(taxon, cluster),
		# records in cluster 
		cluster_n = sum(n),
		# proportion of records in cluster
		cluster_prop = cluster_n/taxon_n,
		# cluster type (major, minor or ND)
		type = dplyr::case_when(
			split == FALSE ~ "major",
			taxon_n < con_min_n ~ "ND",
			# use 'round' to try to avoid floating-point comparison issues when value is either 'exactly' the consensus or its converse
			round(cluster_prop - con_min_prop, 10) >= 0 ~ "major",
			round(cluster_prop - (1 - con_min_prop), 10) <= 0 ~ "minor",
			.default = "ND"
		)
	) %>% 
	# force ND on all clusters for a taxon if no major cluster exists
	dplyr::mutate(
		.by = c(taxon),
		nd = !any(type == "major")
	) %>%
	dplyr::mutate(
		.by = c(taxon, cluster),
		type = dplyr::if_else(
			nd == TRUE, 
			"ND",
			type
		)
	) %>% 
	dplyr::select(-nd) %>%
	# does the genus have consistent internal taxonomy (after filtering)?
	dplyr::mutate(
		.by = c(taxon),
		consistent = "major" %in% unique(type)
	) 
	
# names of sequences that are newly minor-ised (new intragenus outliers)
minor_new <- intragenus_new %>% dplyr::filter(type == "minor", name %in% names(seqs_max)) %>% dplyr::pull(name)

# names of sequences that were removed in max threshold detection or new intragenus filtering
removed_new <- c(minor_new, max_removed)

# recalculate central sequence for all consistent genera
consistent_genera <- 
	intragenus_new %>%
	dplyr::filter(consistent) %>%
	dplyr::pull(taxon) %>%
	unique()

consistent_central <- 
	lapply(
		seq_along(aligned_genera_list),
		function(y){
			x <- aligned_genera_list[[y]]
			# remove any sequences if they have been excluded
			x.ret <- x[!names(x) %in% removed_new]
			# get genus name
			genus_name <- names(x.ret)[1] %>% stringr::str_extract(., "(?<=;).+$") %>% stringr::str_remove(., ";[^;]+$")
			if (genus_name %in% consistent_genera){
				if(length(x.ret) > 1){
					# create distance matrix
					x_d <- DECIPHER::DistanceMatrix(x.ret, verbose = F)
					# median pairwise distance for each sequence
					x_d.median <- apply(x_d, 1, median)
					# indices and names of sequences with smallest median distance
					x_d.min <- which(x_d.median == min(x_d.median))
					# find longest central sequence
					central_idx <- x.ret[names(x.ret) %in% names(x_d.min)] %>% DECIPHER::RemoveGaps(.) %>% lengths() %>% which.max()
					# get name of central sequence
					central_name <- x.ret[names(x.ret) %in% names(x_d.min)][central_idx] %>% names()
				} else {
					central_name <- names(x.ret)
				}
				# create tibble output that links name to genus
				return(tibble::tibble(name = central_name, taxon = genus_name, central = TRUE))
			} else {
				return(NULL)
			}
				
		}
	) %>%
	dplyr::bind_rows()

# join central sequences to new intragenus tibble
intragenus_central <- 
	intragenus_new %>%
	dplyr::left_join(., consistent_central, by = c("name", "taxon")) %>%
	dplyr::mutate(central = dplyr::if_else(is.na(central), FALSE, TRUE))

## make sure every classified consistent genus has a central sequence
# classified consistent genera
cg_check <- intragenus_central %>% dplyr::filter(consistent) %>% dplyr::pull(taxon) %>% unique()
# genera of central sequences
cs_check <- intragenus_central %>% dplyr::filter(central) %>% dplyr::pull(taxon) %>% unique()

if (!identical(cg_check, cs_check)){
	stop("Not all consistent genera have a central sequence")
}

# summarise consistent genera to taxon, central sequence and size of genus in rep-seqs
cg_class <- 
	intragenus_central %>%
	dplyr::rename(genus = taxon) %>%
	dplyr::filter(consistent) %>%
	dplyr::mutate(
		.by = genus, 
		rep_n = sum(n)
	) %>%
	dplyr::filter(central) %>%
	dplyr::select(genus, name, rep_n)

## update synthetic genera to get new central sequences if needed
# remaining sequences in synthetic genera and whether they are central 
synthetic_seqs <- 
	counts_new %>%
	dplyr::mutate(
		genus = stringr::str_extract(name, "(?<=;).+?$") %>% stringr::str_remove(., ";[^;]+?$")
	) %>%
	dplyr::filter(genus %>% stringr::str_detect(., ";Unclassified\\d+$")) %>%
	dplyr::arrange(genus, desc(n)) %>%
	dplyr::filter(!name %in% max_removed) %>%
	dplyr::mutate(
		central = name %in% synthetic_reps
	) 

# which remaining synthetic genera are missing central sequences?
sg_all <- synthetic_seqs %>% dplyr::pull(genus) %>% unique() # all 
sg_cs <- synthetic_seqs %>% dplyr::filter(central) %>% dplyr::pull(genus) %>% unique() # those with central 
sg_cs.missing <- sg_all[!sg_all %in% sg_cs] # those without central 

# get new central sequences for synthetic genera without them by selecting sequence with the highest rep count
sg_cs.new <- 
	synthetic_seqs %>% 
	dplyr::filter(genus %in% sg_cs.missing) %>%
	dplyr::arrange(genus, desc(n)) %>%
	dplyr::slice(.by = genus, 1) %>%
	dplyr::pull(name)

# update synthetic seqs tibble with new central sequences
synthetic_seqs <- 
	synthetic_seqs %>%
	dplyr::mutate(central = dplyr::if_else(name %in% sg_cs.new, TRUE, central))

# summarise synthetic seqs to genera, central seq and total rep count
cg_synth <- 
	synthetic_seqs %>%
	dplyr::mutate(
		.by = genus, 
		rep_n = sum(n)
	) %>%
	dplyr::filter(central) %>%
	dplyr::select(genus, name, rep_n)

# add classified consistent genera and synthetic genera together for a tibble of all consistent genera
cg_all <- dplyr::bind_rows(cg_class, cg_synth)

# write consistent genera to file
readr::write_csv(cg_all, "cg_all.csv")

# write list of consistent genera to file
readr::write_lines(cg_all$genus, "cg_list.txt")
