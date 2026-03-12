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

# value parameters
component_group_size <- as.numeric(component_group_size)
con_min_n <- as.numeric(con_min_n)

allowed_ranks <- c("kingdom","phylum", "class", "order", "family", "genus", "species")
root_ranks <- c("root", allowed_ranks)
root.to.genus <- root_ranks[-8]

# get thresholds as variables
species_min <- thresholds[thresholds$lsr == "species","min"] %>% as.numeric()
genus_min <- thresholds[thresholds$lsr == "genus","min"] %>% as.numeric()
family_min <- thresholds[thresholds$lsr == "family","min"] %>% as.numeric()
order_min <- thresholds[thresholds$lsr == "order","min"] %>% as.numeric()
class_min <- thresholds[thresholds$lsr == "class","min"] %>% as.numeric()
phylum_min <- thresholds[thresholds$lsr == "phylum","min"] %>% as.numeric()
kingdom_min <- thresholds[thresholds$lsr == "kingdom","min"] %>% as.numeric()


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

# make new counts tibble using names(seqs)
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


## make tibble of genus lineage of each sequence record
seqs_genus <- 
	seqs %>%
	names() %>%
	tibble::as_tibble_col(column_name = "name") %>%
    tidyr::separate(col = name, into = c("seqid", "lineage_string"), sep = ";", extra = "merge", remove = F) %>%
    dplyr::mutate(
        genus = stringr::str_extract(lineage_string, "^.+(?=;)")
    ) %>%
	dplyr::select(name, genus)


### loop through comparator alignments and flag genera that fail relevant min threshold (first sequence is the focal record)
flags_cg <- 
	lapply(
		seq_along(comparators_list),
		function(y){
			x <- comparators_list[[y]]
			# name of focal sequence
			focal_name <- names(x)[1]
			
			x_tax <- 
				stringr::str_split(focal_name, ";", n = 8, simplify = F) %>% 
			    unlist() %>% 
				.[2:8] %>%
			    stringr::str_replace_all(., "^Unclassified$", "UNCLASSIFIED")
			
			## NOTE: Don't need to do lowest classification determination as should have been done in previous process
			
			# # which ranks are classified?
			# x_classified <- 
			#     x_tax %>% 
			# 	stringr::str_detect(., "UNCLASSIFIED|Unclassified", negate = T)
			# # get lowest classified rank for focal sequence (ie. max index)
			# x_lcr_i <- which(root_ranks == root_ranks[max(match(root_ranks[x_classified], root_ranks))] )
			
			# distance matrix
			x.d <- DECIPHER::DistanceMatrix(x, verbose = F)
			
			x.c <- 
				x.d[1,] %>%
				tibble::as_tibble(rownames = "comparator_name") %>%
				dplyr::filter(comparator_name != focal_name) %>%
				dplyr::mutate(
					# convert distance to identity
					fid = 1 - value, 
					value = NULL,
					# get taxonomy of comparators
					species = stringr::str_extract(comparator_name, "(?<=;).+?$"),
					genus = stringr::str_extract(species, ".+(?=;)"),
					family = stringr::str_extract(genus, ".+(?=;)"),
			        order = stringr::str_extract(family, ".+(?=;)"),
			        class = stringr::str_extract(order, ".+(?=;)"),
			        phylum = stringr::str_extract(class, ".+(?=;)"),
			        kingdom = stringr::str_extract(phylum, ".+(?=;)"),
					# get LSR of focal and each comparator
					lsr = dplyr::case_when(
				        base::paste(x_tax[1:5], collapse = ";") == family ~ "family",
				        base::paste(x_tax[1:4], collapse = ";") == order ~ "order",
				        base::paste(x_tax[1:3], collapse = ";") == class ~ "class",
				        base::paste(x_tax[1:2], collapse = ";") == phylum ~ "phylum",
				        base::paste(x_tax[1], collapse = ";") == kingdom ~ "kingdom",
				        # shouldn't be these last two so put at the bottom to speed up case matching
				        base::paste(x_tax[1:7], collapse = ";") == species ~ "species",
				        base::paste(x_tax[1:6], collapse = ";") == genus ~ "genus",
				        .default = "root"
					),
					# NOTE: don't need to check lowest classified rank because they shouldn't have been selected as comparators otherwise
					# flag comparator if fid is lower than threshold for LSR
					flag = dplyr::case_when(
						lsr == "family" & fid < family_min ~ TRUE,
						lsr == "order" & fid < order_min ~ TRUE,
						lsr == "class" & fid < class_min ~ TRUE,
						lsr == "phylum" & fid < phylum_min ~ TRUE,
						lsr == "kingdom" & fid < kingdom_min ~ TRUE,
						lsr == "species" & fid < species_min ~ TRUE,
						lsr == "genus" & fid < genus_min ~ TRUE,
						.default = FALSE
					)		
				) 
			
			# there shouldn't be any comparators with species or genus LSR to focal sequence
			if ( any(x.c$lsr %in% c("species","genus"))) stop("One or more comparators have an LSR to focal sequence of species or genus")	
			
			# keep flagged comparators and convert to genus pairs
			x.flags <- 
				x.c %>%
				dplyr::filter(flag) %>%
				dplyr::mutate(
					f_g = base::paste(x_tax[1:6], collapse = ";"),
					focal_name = focal_name
				) 
			
			if (y %% 100 == 0) message(stringr::str_glue("Finished {y} of {length(comparators_list)}"))	
			
			return(x.flags)	
		}
	) %>% 
	dplyr::bind_rows()

## determine which LCGs were not flagged when focal

flagged_focal_cg <- flags_cg$f_g %>% unique()

lcg_nf <- 
	cg_all %>%
	# get LCG	
	dplyr::filter(rep_n >= 2*con_min_n) %>%
	# keep only those that are not flagged as focal
	dplyr::filter(!genus %in% flagged_focal_cg) 

lcg_nf_tax <-
	lcg_nf %>%
	dplyr::mutate(
		genus = paste0("Root;",genus),
		family = stringr::str_extract(genus, ".+(?=;)"),
        order = stringr::str_extract(family, ".+(?=;)"),
        class = stringr::str_extract(order, ".+(?=;)"),
        phylum = stringr::str_extract(class, ".+(?=;)"),
        kingdom = stringr::str_extract(phylum, ".+(?=;)"),
		root = "Root"
	)
	

### for each inconsistent genus, select a single non-flagged LCG (with replacement) per LSR and flag

## determine inconsistent genera from absence in cg_all tibble
# get sequences and genera
ig_ind <- 
	names(seqs) %>%
	tibble::as_tibble_col(column_name = "name") %>%
	dplyr::mutate(
		species = stringr::str_extract(name, "(?<=;).+?$"),
		genus = stringr::str_remove(species, ";[^;]+?$"),
		family = stringr::str_remove(genus, ";[^;]+?$"),
		order = stringr::str_remove(family, ";[^;]+?$"),
		class = stringr::str_remove(order, ";[^;]+?$"),
		phylum = stringr::str_remove(class, ";[^;]+?$"),
		kingdom = stringr::str_remove(phylum, ";[^;]+?$")
	) %>%
	dplyr::filter(!genus %in% cg_all$genus) 

# get just genera
ig_gen <- 
	ig_ind %>%
	dplyr::distinct(genus, .keep_all = TRUE) %>%
	dplyr::pull(genus) 

## for each IG, find a single non-flagged LCG per rank

flags_ig <- 
	lapply(
		seq_along(ig_gen),
		function(y){
			x <- ig_gen[y]
			x_tax <- 
				stringr::str_split(x, ";", n = 7, simplify = F) %>% 
			    unlist() %>% 
				c("Root", .) %>%
			    stringr::str_replace_all(., "^Unclassified$", "UNCLASSIFIED")
			# which ranks are classified?
			x_classified <- 
			    x_tax %>% 
				stringr::str_detect(., "UNCLASSIFIED|Unclassified", negate = T)
			# get lowest classified rank for focal sequence (ie. max index)
			x_lcr_i <- which(root.to.genus == root.to.genus[max(match(root.to.genus[x_classified], root.to.genus))] )
			# select flags
			out <- 
				lcg_nf_tax %>%
				dplyr::mutate(
					lsr = dplyr::case_when(
			            base::paste(x_tax[1:6], collapse = ";") == family ~ "family",
			            base::paste(x_tax[1:5], collapse = ";") == order ~ "order",
			            base::paste(x_tax[1:4], collapse = ";") == class ~ "class",
			            base::paste(x_tax[1:3], collapse = ";") == phylum ~ "phylum",
			            base::paste(x_tax[1:2], collapse = ";") == kingdom ~ "kingdom",
			            .default = "root"
					),
					lcr_i = dplyr::case_when(
			            stringr::str_detect(genus, "Unclassified", negate = T) ~ 7,
			            stringr::str_detect(family, "Unclassified", negate = T) ~ 6,
			            stringr::str_detect(order, "Unclassified", negate = T) ~ 5,
			            stringr::str_detect(class, "Unclassified", negate = T) ~ 4,
			            stringr::str_detect(phylum, "Unclassified", negate = T) ~ 3,
			            stringr::str_detect(kingdom, "Unclassified", negate = T) ~ 2,
			            .default = 1
			        ),
					lcr = pmin(x_lcr_i,lcr_i)
				) %>%
				# subset to genera that are known to be LSR at each rank
				dplyr::filter(
					(lsr == "family" & lcr %in% c(7)) |
					(lsr == "order" & lcr %in% c(7,6)) |
					(lsr == "class" & lcr %in% c(7,6,5)) |
					(lsr == "phylum" & lcr %in% c(7,6,5,4)) |
					(lsr == "kingdom" & lcr %in% c(7,6,5,4,3))
				) %>%
				dplyr::group_by(lsr) %>%
				# select 1 genus per rank
				dplyr::slice_sample(n = 1) %>%
				dplyr::mutate(f_g = x) %>%
				dplyr::ungroup()
			
			return(out)
		}
	) %>% 
	dplyr::bind_rows() %>%
	dplyr::arrange(f_g, lsr) %>%
	dplyr::mutate(
		c_g = stringr::str_remove(genus, "^Root;")
	) 


## convert cg and ig flags tibbles into pairs of genera
# get final flags for CG
flags_cg.f <- 
	flags_cg %>%
	# # select a single flag per LSR per focal genus (with lowest fid)
	# dplyr::arrange(fid) %>%
	# dplyr::group_by(f_g, lsr) %>%
	# dplyr::slice(1) %>%
	# dplyr::ungroup() %>%
	dplyr::select(f_g, c_g = genus)

# get final flags for IG
flags_ig.f <- 
	flags_ig %>%
	dplyr::select(f_g, c_g)

# all flags together
flagged_genera <- dplyr::bind_rows(flags_cg.f, flags_ig.f) %>% dplyr::select(query = f_g, target = c_g)

# build graph of min flagged genera
min_graph <- 
	flagged_genera %>%
	igraph::graph_from_data_frame(
        ., 
        directed = F
    ) %>%
    igraph::simplify() 

message("Built graph")

# determine each component of graph
min_components <- 
    min_graph %>%
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

message("Determined components")

## for components above component_group_size (sum(n_act)), split them into groups (potentially redundant to keep all pairs together)

# split tibble into large and small components
mc_large <- 
	min_components %>%
	dplyr::filter(component_n > component_group_size)

readr::write_csv(mc_large, "mc_large.csv")

mc_small <- 
	min_components %>%
	dplyr::filter(component_n <= component_group_size)

readr::write_csv(mc_small, "mc_small.csv")

message("Defined small and and large components")

# split large component tibble by component, then parse into list format
# ultimately large components, if possible, get split into subcomponents, which are potentially overlapping/redundant connected groups of genera that are m
# some subcomponents will be larger than the max component size if individual genera pairs are by themselves larger than this limit
if (nrow(mc_large) > 1){
	
	# split by component
	mcl_split <- split(mc_large, mc_large$component)
	
	message("Splitting components...")

	# make a tibble of flagged pairs with f_g, c_g, query act_n and target act_n
	flagged_joined <- 
		flagged_genera %>%
		dplyr::left_join(., genus_summary, by = dplyr::join_by(query == genus)) %>%
		dplyr::select(query, target, q_n = act_n) %>%
		dplyr::left_join(., genus_summary, by = dplyr::join_by(target == genus)) %>%
		dplyr::select(query, target, q_n, t_n = act_n) 
		
	## for each component, process, then save all as subcomponents tibble
	mc_sub <- 
		lapply(
			names(mcl_split), # each element in the list is named by the component (eg. '6')
			function(component_id){
				
				message(stringr::str_glue("Component {component_id}, {match(component_id, names(mcl_split))} of {length(mcl_split)}:"))

				# get genera in component
				g_vec <- mcl_split[[component_id]]$genus
				
				# create list of genera pairs, with each genus named with its act_n sequence size
				gp_list <- 
					flagged_joined %>% 
					dplyr::filter(query %in% g_vec, target %in% g_vec) %>%
					dplyr::mutate(
						vec = mapply(c, query, target, q_n, t_n, SIMPLIFY = F)
					) %>%
					dplyr::pull(vec) %>%
					unname() %>% 
					lapply(
						.,
						function(x){
							out <- x[1:2]
							names(out) <- x[3:4]
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
				subcomp_list <- list("0" = NULL)
				for (i in seq_along(gp_list)){
					j <- gp_list[[i]]
					g <- rep(seq_along(subcomp_list), sapply(subcomp_list, length))
					gm <- g[match(j, unlist(subcomp_list))]
					if (all(is.na(gm))){
						# if both genera weren't found, append to list	
						subcomp_list <- base::append(subcomp_list, list(j))
						loc <- length(subcomp_list)
						# make sum of sequences name of new element
						names(subcomp_list)[loc] <- sum(as.numeric(names(j)))
					} else {
						# if one or both genera were found, append the result within the earliest matching element
						first_match <- min(gm[!is.na(gm)])
						updated_match <- c(subcomp_list[[first_match]], j)
						subcomp_list[[first_match]] <- updated_match[!duplicated(updated_match)] # remove repeated genera within element
						loc <- first_match
						# make sum of sequences name of new element
						names(subcomp_list)[loc] <- sum(as.numeric(names(subcomp_list[[first_match]])))
					}
					## use the names of the elements to determine where the newest element should go 
					## this avoids sorting the list each iteration, as you just need to move the focal element so increasing order is maintained
					# find first element larger than the new or updated element
					first_larger <- Position(function(x) x > as.numeric(names(subcomp_list)[loc]), as.numeric(names(subcomp_list)))
					if (is.na(first_larger)){
						# if there was no larger element, move focal element to the end of the list
						subcomp_list <- append(subcomp_list[-loc], subcomp_list[loc], length(subcomp_list))
					} else {
						# if there was a larger element, move focal element to the position before it
						subcomp_list <- append(subcomp_list[-loc], subcomp_list[loc], first_larger)
					}
					
					#subcomp_list <- subcomp_list[order(seq_sum, method = "radix")]
					if ((i %% 100 == 0) || (i == length(gp_list))) message(stringr::str_glue("\t{i} of {length(gp_list)}"))
				}

				# remove remaining null elements
				subcomp_list[sapply(subcomp_list, is.null)] <- NULL
				
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
		component = character(),
		act_n = integer(),
		component_n = integer()
	)
}

# combine small and subcomponented-large components together
mc_all <- dplyr::bind_rows(mc_small, mc_sub)

message("Combined small and processed large components")

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

message("Defined component groups")

# split by component group
mc_split <- 
	mc_all %>%
	dplyr::left_join(., component_groups, by = "component") %>%
	dplyr::arrange(group, component) %>% 
	split(., .$group) 

readr::write_csv(mc_split %>% dplyr::bind_rows(), "component_groups.csv")

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

#stop()