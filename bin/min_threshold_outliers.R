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
	"seqs_file",
    "rf_counts_tsv",
    "thresholds_csv",
    "con_min_n",
    "con_min_prop"
)
lapply(nf_vars, nf_var_check)

### process variables 

# read in list of sequences
alignment_single <- readr::read_lines(alignment_cfa)

# process concatenated fasta format into list of DSS format
alignment_list <- 
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

# read in all sequences
seqs <- Biostrings::readDNAStringSet(seqs_file)

# read in rf_counts csv
counts_new <- readr::read_tsv(rf_counts_tsv, show_col_types = FALSE)

# read in thresholds csv
thresholds <- readr::read_csv(thresholds_csv,  show_col_types = FALSE)

# consensus parameters
con_min_n <- as.numeric(con_min_n)
con_min_prop <- as.numeric(con_min_prop)

### run code

## convert min thresholds to distances
t_i <- thresholds %>%
	dplyr::mutate(threshold_name = paste0(lsr,"_min")) %>%
	dplyr::select(threshold_name, min) %>%
	tibble::deframe()

t_d <- 1 - t_i 

## loop through component groups and identity minor sequences
min_loop <-
	lapply(
        seq_along(alignment_list),
        function(y){
			# get distance matrix from alignment
			x <- alignment_list[[y]]
			x.d <- DECIPHER::DistanceMatrix(x, verbose = F)
			
			# loop through the thresholds, producing clusters for each one and then combining into a single tibble
			min_cl <- 
			    lapply(
			    	names(t_d),
			        function(y){
			        	DECIPHER::TreeLine(
							myDistMatrix = x.d,
			                method = "single",
			                type = "clusters",
			                cutoff = t_d[names(t_d) == y],
			                verbose = F
			            ) %>%
			                tibble::as_tibble(rownames = "name") %>%
			                tidyr::pivot_longer(cluster, names_to = "threshold", values_to = "cluster") %>%
			                dplyr::mutate(threshold = replace(threshold, threshold == "cluster", y))
			        }
			    ) %>%
			    dplyr::bind_rows() %>%
			    tidyr::pivot_wider(
			        names_from = "threshold",
			        values_from = "cluster"
			    ) %>%
			    # get taxonomy at each rank for each cluster
			    dplyr::mutate(
			        # add component group hash
			        component_group = rlang::hash(names(x)),
			        species = stringr::str_extract(name, "(?<=;).+$"),
			        genus = stringr::str_extract(species, ".+(?=;)"),
			        family = stringr::str_extract(genus, ".+(?=;)"),
			        order = stringr::str_extract(family, ".+(?=;)"),
			        class = stringr::str_extract(order, ".+(?=;)"),
			        phylum = stringr::str_extract(class, ".+(?=;)"),
			        kingdom = stringr::str_extract(phylum, ".+(?=;)")
			    ) %>%
			    # get record n for each sequence
			    dplyr::left_join(., counts_new, by = "name") %>%
			    dplyr::arrange(kingdom_min, phylum_min, class_min, order_min, family_min, genus_min, species_min, desc(n))
			
			#### determine mixed clusters
			
			## kingdom
			c_k_min <- 
			    min_cl %>%
			    dplyr::select(component_group, name, n, taxon = kingdom, cluster = kingdom_min) %>%
			    # remove sequences without kingdom classification
			    dplyr::filter(!stringr::str_detect(taxon, ";Unclassified\\d+$")) %>%
			    dplyr::mutate(
					.by = c(taxon),
					threshold = "kingdom_min",
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
			    dplyr::arrange(taxon, desc(cluster_prop), desc(n)) 
			
			# get names of records that failed step 
			k_min_fail <- c_k_min %>% dplyr::filter(type == "minor") %>% pull(name)
			# get names of records that passed step 
			k_min_pass <- c_k_min %>% dplyr::filter(type %in% c("major", "ND")) %>% pull(name)
			
			
			## phylum
			c_p_min <- 
			    min_cl %>%
			    dplyr::select(component_group, name, n, taxon = phylum, cluster = phylum_min) %>%
			    # remove sequences without phylum classification
			    dplyr::filter(!stringr::str_detect(taxon, ";Unclassified\\d+$")) %>%
				# remove sequences previously removed
			    dplyr::filter(!name %in% c(k_min_fail)) %>%
			    dplyr::mutate(
					.by = c(taxon),
					threshold = "phylum_min",
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
			    dplyr::arrange(taxon, desc(cluster_prop), desc(n)) 
			
			# get names of records that failed step 
			p_min_fail <- c_p_min %>% dplyr::filter(type == "minor") %>% pull(name)
			# get names of records that passed step 
			p_min_pass <- c_p_min %>% dplyr::filter(type %in% c("major", "ND")) %>% pull(name)
			
			
			## class
			c_c_min <- 
			    min_cl %>%
			    dplyr::select(component_group, name, n, taxon = class, cluster = class_min) %>%
			    # remove sequences without class classification
			    dplyr::filter(!stringr::str_detect(taxon, ";Unclassified\\d+$")) %>%
				# remove sequences previously removed
			    dplyr::filter(!name %in% c(k_min_fail, p_min_fail)) %>%
			    dplyr::mutate(
					.by = c(taxon),
					threshold = "class_min",
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
			    dplyr::arrange(taxon, desc(cluster_prop), desc(n)) 
			
			# get names of records that failed step 
			c_min_fail <- c_c_min %>% dplyr::filter(type == "minor") %>% pull(name)
			# get names of records that passed step 
			c_min_pass <- c_c_min %>% dplyr::filter(type %in% c("major", "ND")) %>% pull(name)
			
			
			## order
			c_o_min <- 
			    min_cl %>%
			    dplyr::select(component_group, name, n, taxon = order, cluster = order_min) %>%
			    # remove sequences without order classification
			    dplyr::filter(!stringr::str_detect(taxon, ";Unclassified\\d+$")) %>%
				# remove sequences previously removed
			    dplyr::filter(!name %in% c(k_min_fail, p_min_fail, c_min_fail)) %>%
			    dplyr::mutate(
					.by = c(taxon),
					threshold = "order_min",
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
			    dplyr::arrange(taxon, desc(cluster_prop), desc(n)) 
			
			# get names of records that failed step 
			o_min_fail <- c_o_min %>% dplyr::filter(type == "minor") %>% pull(name)
			# get names of records that passed step 
			o_min_pass <- c_o_min %>% dplyr::filter(type %in% c("major", "ND")) %>% pull(name)
			
			
			## family
			c_f_min <- 
			    min_cl %>%
			    dplyr::select(component_group, name, n, taxon = family, cluster = family_min) %>%
			    # remove sequences without family classification
			    dplyr::filter(!stringr::str_detect(taxon, ";Unclassified\\d+$")) %>%
				# remove sequences previously removed
			    dplyr::filter(!name %in% c(k_min_fail, p_min_fail, c_min_fail, o_min_fail)) %>%
			    dplyr::mutate(
					.by = c(taxon),
					threshold = "family_min",
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
			    dplyr::arrange(taxon, desc(cluster_prop), desc(n)) 
			
			# get names of records that failed step 
			f_min_fail <- c_f_min %>% dplyr::filter(type == "minor") %>% pull(name)
			# get names of records that passed step 
			f_min_pass <- c_f_min %>% dplyr::filter(type %in% c("major", "ND")) %>% pull(name)
			
			
			
			## combined output for component_group
			c_all_min <-
				dplyr::bind_rows(
					c_k_min,
					c_p_min,
					c_c_min,
					c_o_min,
					c_f_min
				)
			
			## get tibble of minor sequences and their violated thresholds
			combined_minor <- 
				c_all_min %>% 
				dplyr::filter(type == "minor") %>%
				dplyr::select(component_group, name, n, threshold)
			
			message(paste("Finished component",y))
			
			return(
				list(
					"all" = c_all_min,
					"minor" = combined_minor
				)
			)
        }
	)

min_all <- lapply(min_loop, '[[', 1) %>% dplyr::bind_rows()
min_minor <- lapply(min_loop, '[[', 2) %>% dplyr::bind_rows()

names_minor <- min_minor$name %>% unique()

seqs_retained <- seqs[!names(seqs) %in% names_minor]
seqs_removed <- seqs[names(seqs) %in% names_minor]

# output retained and removed sequences
if (length(seqs_retained) > 0){
	Biostrings::writeXStringSet(seqs_retained, "retained.fasta", width = 9999)
} else {
	file.create("retained.fasta")
}

if (length(seqs_removed) > 0){
	Biostrings::writeXStringSet(seqs_removed, "removed.fasta", width = 9999)
} else {
	file.create("removed.fasta")
}

# output tables of filtering info
readr::write_csv(min_all, "min_tibble.csv")
readr::write_csv(min_minor, "min_minor.csv")