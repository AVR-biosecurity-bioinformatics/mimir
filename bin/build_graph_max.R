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

flagged_genera <- readr::read_csv(flagged_genera_file, show_col_types = FALSE) %>% dplyr::distinct()

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

# build graph of max violations (connected genera)
max_graph <- 
	flagged_genera %>%
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
	# get size of each genus
	dplyr::left_join(., genus_summary, by = "genus") %>%
	dplyr::select(-rep_n) %>%
	# get size of each component
	dplyr::mutate(
		.by = component,
		component_n = sum(act_n),
		component = as.character(component) # this is needed later combining with subcomponent ids (of form "x_y")
	) 
	
# for components above component_group_size (sum(n_act)), split them into groups (potentially redundant to keep all pairs together)

# split tibble into large and small components
mc_large <- 
	max_components %>%
	dplyr::mutate(
		.by = component,
		component_n = sum(act_n)
	) %>%
	dplyr::filter(component_n > component_group_size)

mc_small <- 
	max_components %>%
	dplyr::mutate(
		.by = component,
		component_n = sum(act_n)
	) %>%
	dplyr::filter(component_n <= component_group_size)
	
# split large component tibble by component, then parse into list format
# ultimately large components, if possible, get split into subcomponents, which are potentially overlapping/redundant connected groups of genera that are m
# some subcomponents will be larger than the max component size if individual genera pairs are by themselves larger than this limit
if (nrow(mc_large) > 1){
	
	# split by component
	mcl_split <- split(mc_large, mc_large$component)
	
	## for each component, process, then save all as subcomponents tibble
	mc_sub <- 
		lapply(
			names(mcl_split), # each element in the list is named by the component (eg. '6')
			function(component_id){
				
				# get genera in component
				g_vec <- mcl_split[[component_id]]$genus
				
				# create list of genera pairs, with each genus named with its act_n sequence size
				gp_list <- 
					flagged_genera %>% 
					dplyr::filter(query %in% g_vec, target %in% g_vec) %>%
					dplyr::mutate(
						vec = mapply(c, query, target, SIMPLIFY = F)
					) %>%
					dplyr::pull(vec) %>%
					unname() %>%
					lapply(
						.,
						function(x){
							out <- x
							names(out)[1] <- genus_summary[genus_summary$genus == x[1],]$act_n
							names(out)[2] <- genus_summary[genus_summary$genus == x[2],]$act_n
							return(out)
						}
					) 

				# reorder list of pairs so the pairs with smallest total sequences are at the start
				gp_order <- 
					gp_list %>%
					lapply(
						.,
						function(x){
							sum(as.numeric(names(x)))
						}
					) %>%
					unlist() %>%
					order()
				
				gp_list <- gp_list[gp_order]
				
				# group pairs into potentially redundant subcomponents
				subcomp_list <- list()
				for (i in seq_along(gp_list)){
					found <- FALSE
					# append first instance of matching genus to that element of subcomp_list
					for (j in seq_along(subcomp_list)){
						if (any(gp_list[[i]] %in% subcomp_list[[j]])){
							found <- TRUE
							new_j <- base::append(subcomp_list[[j]], gp_list[[i]])
							subcomp_list[[j]] <- new_j[!duplicated(new_j)] # remove repeated genera within element
							break()
						}
					}
					# if neither genus was found in any existing elements, create a new element in subcomp_list
					if (!found){
						subcomp_list <- base::append(subcomp_list, gp_list[i])
					}
					# arrange subcomp_list in ascending order of total number of sequences in the genera, so smallest groups are compared first
					seq_sum <- sapply(subcomp_list, function(x) names(x) %>% as.numeric() %>% sum())
					subcomp_list <- subcomp_list[order(seq_sum, decreasing = F)]
				}

				# convert subcomp_list into a tibble format (implicit output)
				lapply(
					seq_along(subcomp_list), 
					function(idx){
						x <- subcomp_list[[idx]]
						out <- 
							x %>% 
							tibble::enframe(., name = "act_n", value = "genus") %>%
							dplyr::mutate(
								act_n = as.numeric(act_n),
								component_n = sum(act_n),
								# make new component name using subcomponent id and original component id
								component = paste0(component_id, "_", idx)
							) %>%
							dplyr::select(genus, component, act_n, component_n)
						return(out)
					}
				) %>%
					dplyr::bind_rows()
			}
		) %>%
		dplyr::bind_rows()
			
} else {
	# if there were no large components to split, output empty subcomponents tibble
	mc_sub <- tibble::tibble(
		genus = character(),
		component = double(),
		act_n = integer(),
		component_n = integer()
	)
}

# combine small and subcomponented-large components together
mc_all <- dplyr::bind_rows(mc_small, mc_sub)

# group small components together so groups with more than 1 member have max 'component_group_size' sequences
# this is so we don't align more than 'component_group_size' sequences together unless we have to
component_groups <- 
	mc_all %>%
	dplyr::distinct(component, component_n) %>%
	dplyr::arrange(component_n) %>%
	dplyr::mutate(
		# idea from: https://stackoverflow.com/questions/34531568/conditional-cumsum-with-reset
		cs = base::Reduce(\(x, y) if (x + y <= component_group_size) x + y else y, x = component_n, accumulate = TRUE),
		group = cumsum(component_n == cs)
	) %>%
	dplyr::select(component, group) 

# split by component group
mc_split <- 
	mc_all %>%
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
			dplyr::filter(genus %in% x$genus) %>% # repeated genera names don't affect this selection method
			dplyr::pull(name) 
		x_seqs <- seqs[names(seqs) %in% x_seq_names]
		Biostrings::writeXStringSet(x_seqs, stringr::str_glue("component_group{y}.fasta"), width = 9999)
		message(stringr::str_glue("{length(x_seqs)} sequences in component group {y}"))
	}
)