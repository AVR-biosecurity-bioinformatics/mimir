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
    "alignment_cfa",
    "rf_counts_tsv",
    "thresholds_csv",
    "con_min_n",
    "con_min_prop"
)
lapply(nf_vars, nf_var_check)

### process variables 

# read in list of sequences
alignment_single <- readr::read_lines(alignment_cfa)

# process concatenated fasta format into DSS format
seqs_list <- 
	alignment_single %>%
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

# read in rf_counts csv
rf_counts <- readr::read_tsv(rf_counts_tsv, col_names = c("name","n"), show_col_types = FALSE)

# read in thresholds csv
thresholds <- readr::read_csv(thresholds_csv,  show_col_types = FALSE)

# get species and genus min thresholds as distances
species_min_d <- 1 - as.numeric(thresholds[thresholds$lsr == "species", "min"])
genus_min_d <- 1 - as.numeric(thresholds[thresholds$lsr == "genus", "min"])

# consensus parameters
con_min_n <- as.numeric(con_min_n)
con_min_prop <- as.numeric(con_min_prop)

### run code

# loop through each DSS object (representing each genus), cluster and detect outliers (as well as defining central sequences)
gs_outliers <- 
	lapply(
		seq_along(seqs_list),
		function(y){
			#browser()
			message(stringr::str_glue("Starting {y} of {length(seqs_list)}"))
			x <- seqs_list[[y]]
			if (length(x) > 1){
				x_d <- DECIPHER::DistanceMatrix(x, verbose = F)
				# produce distance matrix, cluster, then organise into tibble
				g_c <- 
					x_d %>% 
					DECIPHER::TreeLine(
						myDistMatrix = ., 
						method = "single",
						type = "clusters",
						cutoff = c(species_min_d, genus_min_d),
						verbose = F
					) %>%
					tibble::as_tibble(rownames = "name") %>%
					dplyr::rename(species_min = 2, genus_min = 3) %>%
					# add rep counts
					dplyr::left_join(., rf_counts, by = "name") %>%
					# get species name per sequence
					dplyr::mutate(
						species = stringr::str_extract(name, "(?<=;).+?$"),
						genus = stringr::str_remove(species, ";[^;]+?$")
					) %>%
					dplyr::arrange(genus_min, species_min, desc(n))
			} else {
				# create dummy cluster output when there is only a single sequence in the genus
				g_c <- 
					tibble(name = names(x), species_min = 1, genus_min = 1) %>%
					dplyr::left_join(., rf_counts, by = "name") %>%
					# get species name per sequence
					dplyr::mutate(
						species = stringr::str_extract(name, "(?<=;).+?$"),
						genus = stringr::str_remove(species, ";[^;]+?$")
					) %>%
					dplyr::arrange(genus_min, species_min, desc(n))
			}


			## determine if genus is split into multiple clusters, and whether they are major, minor or ND
			## also determine if genus is consistent 
			genus_check <- 
				g_c %>%
				# remove sequences without genus classification
				dplyr::filter(!stringr::str_detect(genus, ";Unclassified$")) %>%
				dplyr::select(name, n, taxon = genus, cluster = genus_min) %>%
				dplyr::group_by(taxon) %>%
				dplyr::mutate(
					threshold = "genus_min",
					# is taxon split?
					split = dplyr::if_else(length(unique(cluster)) == 1, FALSE, TRUE),
					# records in genus
					taxon_n = sum(n)
				) %>%
				dplyr::group_by(taxon, cluster) %>%
				dplyr::mutate(
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
				dplyr::group_by(taxon) %>%
				dplyr::mutate(
					nd = !any(type == "major")
				) %>%
				dplyr::group_by(taxon, cluster) %>%
				dplyr::mutate(
					type = dplyr::if_else(
						nd == TRUE, 
						"ND",
						type
					)
				) %>%
				dplyr::select(-nd) %>%
				# does the genus have consistent internal taxonomy (after filtering)?
				dplyr::group_by(taxon) %>%
				dplyr::mutate(consistent = "major" %in% unique(type)) %>%
				dplyr::ungroup()
			
			# get DSS of sequences removed at genus level
			genus_minor_seqs <- x[names(x) %in% (genus_check %>% dplyr::filter(type == "minor") %>% .$name)]
			
			# check for split species with remaining sequences
			species_check <- 
				g_c %>%
					# remove sequences without species classification
				dplyr::filter(!stringr::str_detect(species, ";Unclassified$")) %>%
					# remove sequences removed at genus level
					dplyr::filter(!name %in% names(genus_minor_seqs)) %>%
							dplyr::select(name, n, taxon = species, cluster = species_min) %>%
					dplyr::arrange(taxon, cluster, desc(n)) %>%
				dplyr::group_by(taxon) %>%
				dplyr::mutate(
						threshold = "species_min",
					# is taxon split?
					split = dplyr::if_else(length(unique(cluster)) == 1, FALSE, TRUE),
					# records in species
					taxon_n = sum(n)
				) %>%
				dplyr::group_by(taxon, cluster) %>%
				dplyr::mutate(
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
				# force ND on all clusters for a species if no major cluster exists
				dplyr::group_by(taxon) %>%
				dplyr::mutate(
					nd = !any(type == "major")
				) %>%
				dplyr::group_by(taxon, cluster) %>%
				dplyr::mutate(
					type = dplyr::if_else(
						nd == TRUE, 
						"ND",
						type
					)
				) %>%
				dplyr::select(-nd) %>%
				# does the species have consistent internal taxonomy (after filtering)?
				dplyr::group_by(taxon) %>%
				dplyr::mutate(consistent = "major" %in% unique(type)) %>%
				dplyr::ungroup() 
			
			# get DSS of sequences removed at species level
			species_minor_seqs <- x[names(x) %in% (species_check %>% dplyr::filter(type == "minor") %>% .$name)]
			
			# names of all removed sequences
			minor_names <- c(names(species_minor_seqs), names(genus_minor_seqs))
			
			# combine tibbles
			out_tibble.pre <- 
				dplyr::bind_rows(genus_check, species_check)
			
			# determine central sequence of genus (if consistent)
			if (all(genus_check$consistent)){
				# only use distance matrix if there is more than one sequence
				if(length(x) > 1){
					# remove minor sequences from the distance matrix
					x_d.new <- x_d[!rownames(x_d) %in% minor_names, !colnames(x_d) %in% minor_names, drop = F]
					if(length(rownames(x_d.new)) == 1){
						# if only a single sequence remains, use that as central sequence
						central_name <- rownames(x_d.new)
					} else {
						# median pairwise distance for each sequence
						x_d.median <- apply(x_d.new, 1, median)
						# indices and names of sequences with smallest median distance
						x_d.min <- which(x_d.median == min(x_d.median))
						# find longest central sequence
						central_idx <- x[names(x) %in% names(x_d.min)] %>% DECIPHER::RemoveGaps(.) %>% lengths() %>% which.max()
						# get name of central sequence
						central_name <- x[names(x) %in% names(x_d.min)][central_idx] %>% names()
					}
				} else {
					central_name <- names(x)
				}
				
				# add central flag to out_tibble
				out_tibble <- 
					out_tibble.pre %>%
					dplyr::mutate(central = name == central_name)
				 
			} else {
				# don't define central sequence
				out_tibble <- 
					out_tibble.pre %>%
					dplyr::mutate(central = FALSE)
			}
			
			### return list per genus
			return(
				list(
					"tibble" = out_tibble,
					"genus_minor_seqs" = genus_minor_seqs,
					"species_minor_seqs" = species_minor_seqs
				)
			)
	
		}
	)


gs_tibble <- lapply(gs_outliers, '[[', 1) %>% dplyr::bind_rows()
gs_gminor <- lapply(gs_outliers, '[[', 2) %>% do.call(c, .)
gs_sminor <- lapply(gs_outliers, '[[', 3) %>% do.call(c, .)

readr::write_csv(gs_tibble, "gs_tibble.csv")

if (length(gs_gminor) > 0){
	Biostrings::writeXStringSet(gs_gminor, "gminor.fasta", width = 9999)
} else {
	file.create("gminor.fasta")
}

if (length(gs_sminor) > 0){
	Biostrings::writeXStringSet(gs_sminor, "sminor.fasta", width = 9999)
} else {
	file.create("sminor.fasta")
}

seqs_retained <- 
	seqs_list %>% 
	do.call(c, .) %>%
	DECIPHER::RemoveGaps(., removeGaps = "all") %>%
	.[!names(.) %in% c(names(gs_gminor), names(gs_sminor))]

if (length(seqs_retained) > 0){
	Biostrings::writeXStringSet(seqs_retained, "retained.fasta", width = 9999)
} else {
	file.create("retained.fasta")
}