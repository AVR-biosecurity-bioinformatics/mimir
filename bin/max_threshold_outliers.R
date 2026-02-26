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
rf_counts <- readr::read_tsv(rf_counts_tsv, col_names = c("name","n"), show_col_types = FALSE)

# read in thresholds csv
thresholds <- readr::read_csv(thresholds_csv,  show_col_types = FALSE)

# consensus parameters
con_min_n <- as.numeric(con_min_n)
con_min_prop <- as.numeric(con_min_prop)

### run code

## make sure counts sequence names match new names for synthetic genera
# LCR family or above
counts_lcrf <- 
	rf_counts %>%
	dplyr::filter(stringr::str_detect(name, ";Unclassified;[^;]+$"))

# LCR genus or below
count_lcrg <- 
	rf_counts %>%
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

## convert max thresholds to distances
t_i <- thresholds %>%
	dplyr::mutate(threshold_name = paste0(lsr,"_max")) %>%
	dplyr::select(threshold_name, max) %>%
	tibble::deframe()

t_d <- 1 - t_i 

## loop through component groups and identity minor sequences
max_loop <-
	lapply(
        seq_along(alignment_list),
        function(y){
			# get distance matrix from alignment
			x <- alignment_list[[y]]
			
			# loop through the thresholds, producing clusters for each one and then combining into a single tibble
			max_cl <- 
			    lapply(
			    	names(t_d),
			        function(y){
						DECIPHER::DistanceMatrix(x, verbose = F) %>%
			        	DECIPHER::TreeLine(
							myDistMatrix = .,
			                method = "complete",
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
			    dplyr::arrange(kingdom_max, phylum_max, class_max, order_max, family_max, genus_max, species_max, desc(n))
			
			
			#### determine mixed clusters
			
			## family
			c_f_max <- 
			    max_cl %>%
			    dplyr::select(component_group, name, n, taxon = genus, cluster = family_max) %>%
			    # remove sequences without genus classification
			    dplyr::filter(!stringr::str_detect(taxon, ";Unclassified\\d+$")) %>%
			    dplyr::group_by(cluster) %>%
			    dplyr::mutate(
			        threshold = "family_max",
			        # is cluster mixed?
			        mixed = dplyr::if_else(length(unique(taxon)) == 1, FALSE, TRUE) ,
			        # records in cluster
			        cluster_n = sum(n)
			    ) %>%
			    group_by(taxon, cluster) %>%
			    dplyr::mutate(
			        # records in taxon 
			        taxon_n = sum(n),
			        # proportion of records in taxon
			        taxon_prop = taxon_n/cluster_n,
			        # taxon type (major, minor or ND)
			        type = dplyr::case_when(
			            mixed == FALSE ~ "major",
			            cluster_n < con_min_n ~ "ND",
			            # use 'round' to try to avoid floating-point comparison issues when value is either 'exactly' the consensus or its converse
			            round(taxon_prop - con_min_prop, 10) >= 0 ~ "major",
			            round(taxon_prop - (1 - con_min_prop), 10) <= 0 ~ "minor",
			            .default = "ND"
			        )
			    ) %>%
			    # force ND on all taxa for a cluster if no major taxon exists
			    dplyr::group_by(cluster) %>%
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
			    dplyr::ungroup() %>%
			    dplyr::arrange(cluster, desc(taxon_prop), desc(n)) 
			
			# get names of records that failed step 
			f_max_fail <- c_f_max %>% dplyr::filter(type == "minor") %>% pull(name)
			# get names of records that passed step 
			f_max_pass <- c_f_max %>% dplyr::filter(type %in% c("major", "ND")) %>% pull(name)
			            
			
			## order
			c_o_max <- 
			    max_cl %>%
			    # retain only previously passing sequences
			    dplyr::filter(!name %in% c(f_max_fail)) %>%
			    # remove sequences without family classification
			    dplyr::filter(!stringr::str_detect(family, ";Unclassified$")) %>%
			    dplyr::select(component_group, name, n, taxon = family, cluster = order_max) %>%
			    dplyr::group_by(cluster) %>%
			    dplyr::mutate(
			        threshold = "order_max",
			        # is cluster mixed?
			        mixed = dplyr::if_else(length(unique(taxon)) == 1, FALSE, TRUE) ,
			        # records in cluster
			        cluster_n = sum(n)
			    ) %>%
			    group_by(taxon, cluster) %>%
			    dplyr::mutate(
			        # records in taxon 
			        taxon_n = sum(n),
			        # proportion of records in taxon
			        taxon_prop = taxon_n/cluster_n,
			        # taxon type (major, minor or ND)
			        type = dplyr::case_when(
			            mixed == FALSE ~ "major",
			            cluster_n < con_min_n ~ "ND",
			            # use 'round' to try to avoid floating-point comparison issues when value is either 'exactly' the consensus or its converse
			            round(taxon_prop - con_min_prop, 10) >= 0 ~ "major",
			            round(taxon_prop - (1 - con_min_prop), 10) <= 0 ~ "minor",
			            .default = "ND"
			        )
			    ) %>%
			    # force ND on all taxa for a cluster if no major taxon exists
			    dplyr::group_by(cluster) %>%
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
			    dplyr::ungroup() %>%
			    dplyr::arrange(cluster, desc(taxon_prop), desc(n)) 
			
			# get names of records that failed step 
			o_max_fail <- c_o_max %>% dplyr::filter(type == "minor") %>% pull(name)
			# get names of records that passed step 
			o_max_pass <- c_o_max %>% dplyr::filter(type %in% c("major", "ND")) %>% pull(name)
			
			## class
			c_c_max <- 
			    max_cl %>%
			    # retain only previously passing sequences
			    dplyr::filter(!name %in% c(f_max_fail, o_max_fail)) %>%
			    # remove sequences without order classification
			    dplyr::filter(!stringr::str_detect(order, ";Unclassified$")) %>%
			    dplyr::select(component_group, name, n, taxon = order, cluster = class_max) %>%
			    dplyr::group_by(cluster) %>%
			    dplyr::mutate(
			        threshold = "class_max",
			        # is cluster mixed?
			        mixed = dplyr::if_else(length(unique(taxon)) == 1, FALSE, TRUE) ,
			        # records in cluster
			        cluster_n = sum(n)
			    ) %>%
			    group_by(taxon, cluster) %>%
			    dplyr::mutate(
			        # records in taxon 
			        taxon_n = sum(n),
			        # proportion of records in taxon
			        taxon_prop = taxon_n/cluster_n,
			        # taxon type (major, minor or ND)
			        type = dplyr::case_when(
			            mixed == FALSE ~ "major",
			            cluster_n < con_min_n ~ "ND",
			            # use 'round' to try to avoid floating-point comparison issues when value is either 'exactly' the consensus or its converse
			            round(taxon_prop - con_min_prop, 10) >= 0 ~ "major",
			            round(taxon_prop - (1 - con_min_prop), 10) <= 0 ~ "minor",
			            .default = "ND"
			        )
			    ) %>%
			    # force ND on all taxa for a cluster if no major taxon exists
			    dplyr::group_by(cluster) %>%
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
			    dplyr::ungroup() %>%
			    dplyr::arrange(cluster, desc(taxon_prop), desc(n)) 
			
			# get names of records that failed step 
			c_max_fail <- c_c_max %>% dplyr::filter(type == "minor") %>% pull(name)
			# get names of records that passed step 
			c_max_pass <- c_c_max %>% dplyr::filter(type %in% c("major", "ND")) %>% pull(name)
			
			## phylum
			c_p_max <- 
			    max_cl %>%
			    # retain only previously passing sequences
			    dplyr::filter(!name %in% c(f_max_fail, o_max_fail, c_max_fail)) %>%
			    # remove sequences without class classification
			    dplyr::filter(!stringr::str_detect(class, ";Unclassified$")) %>%
			    dplyr::select(component_group, name, n, taxon = class, cluster = phylum_max) %>%
			    dplyr::group_by(cluster) %>%
			    dplyr::mutate(
			        threshold = "phylum_max",
			        # is cluster mixed?
			        mixed = dplyr::if_else(length(unique(taxon)) == 1, FALSE, TRUE) ,
			        # records in cluster
			        cluster_n = sum(n)
			    ) %>%
			    group_by(taxon, cluster) %>%
			    dplyr::mutate(
			        # records in taxon 
			        taxon_n = sum(n),
			        # proportion of records in taxon
			        taxon_prop = taxon_n/cluster_n,
			        # taxon type (major, minor or ND)
			        type = dplyr::case_when(
			            mixed == FALSE ~ "major",
			            cluster_n < con_min_n ~ "ND",
			            # use 'round' to try to avoid floating-point comparison issues when value is either 'exactly' the consensus or its converse
			            round(taxon_prop - con_min_prop, 10) >= 0 ~ "major",
			            round(taxon_prop - (1 - con_min_prop), 10) <= 0 ~ "minor",
			            .default = "ND"
			        )
			    ) %>%
			    # force ND on all taxa for a cluster if no major taxon exists
			    dplyr::group_by(cluster) %>%
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
			    dplyr::ungroup() %>%
			    dplyr::arrange(cluster, desc(taxon_prop), desc(n)) 
			
			# get names of records that failed step 
			p_max_fail <- c_p_max %>% dplyr::filter(type == "minor") %>% pull(name)
			# get names of records that passed step 
			p_max_pass <- c_p_max %>% dplyr::filter(type %in% c("major", "ND")) %>% pull(name)
			
			## kingdom
			c_k_max <- 
			    max_cl %>%
			    # retain only previously passing sequences
			    dplyr::filter(!name %in% c(f_max_fail, o_max_fail, c_max_fail, p_max_fail)) %>%
			    # remove sequences without phylum classification
			    dplyr::filter(!stringr::str_detect(phylum, ";Unclassified$")) %>%
			    dplyr::select(component_group, name, n, taxon = phylum, cluster = kingdom_max) %>%
			    dplyr::group_by(cluster) %>%
			    dplyr::mutate(
			        threshold = "kingdom_max",
			        # is cluster mixed?
			        mixed = dplyr::if_else(length(unique(taxon)) == 1, FALSE, TRUE) ,
			        # records in cluster
			        cluster_n = sum(n)
			    ) %>%
			    group_by(taxon, cluster) %>%
			    dplyr::mutate(
			        # records in taxon 
			        taxon_n = sum(n),
			        # proportion of records in taxon
			        taxon_prop = taxon_n/cluster_n,
			        # taxon type (major, minor or ND)
			        type = dplyr::case_when(
			            mixed == FALSE ~ "major",
			            cluster_n < con_min_n ~ "ND",
			            # use 'round' to try to avoid floating-point comparison issues when value is either 'exactly' the consensus or its converse
			            round(taxon_prop - con_min_prop, 10) >= 0 ~ "major",
			            round(taxon_prop - (1 - con_min_prop), 10) <= 0 ~ "minor",
			            .default = "ND"
			        )
			    ) %>%
			    # force ND on all taxa for a cluster if no major taxon exists
			    dplyr::group_by(cluster) %>%
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
			    dplyr::ungroup() %>%
			    dplyr::arrange(cluster, desc(taxon_prop), desc(n)) 
			
			# get names of records that failed step 
			k_max_fail <- c_k_max %>% dplyr::filter(type == "minor") %>% pull(name)
			# get names of records that passed step 
			k_max_pass <- c_k_max %>% dplyr::filter(type %in% c("major", "ND")) %>% pull(name)
			
			## combined output for component_group
			c_all_max <-
				dplyr::bind_rows(
					c_f_max,
					c_o_max,
					c_c_max,
					c_p_max,
					c_k_max
				)
			
			## get tibble of minor sequences and their violated thresholds
			combined_minor <- 
				c_all_max %>% 
				dplyr::filter(type == "minor") %>%
				dplyr::select(component_group, name, n, threshold)
			
			message(paste("Finished component",y))
			
			return(
				list(
					"all" = c_all_max,
					"minor" = combined_minor
				)
			)
        }
	)

max_all <- lapply(max_loop, '[[', 1) %>% dplyr::bind_rows()
max_minor <- lapply(max_loop, '[[', 2) %>% dplyr::bind_rows()

names_minor <- max_minor$name %>% unique()

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
readr::write_csv(max_all, "max_tibble.csv")
readr::write_csv(max_minor, "max_minor.csv")

