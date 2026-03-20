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
    "cg_csv",
	"cg_list",
    "seqs_max_fasta",
    "con_min_n",
    "con_min_prop"
)
lapply(nf_vars, nf_var_check)

### process variables 

# read in intragenus csv
cg_all <- readr::read_csv(cg_csv, show_col_types = FALSE)

# read in vector of consistent genera to find comparators for (add Root too)
cg_vec <- readr::read_lines(cg_list) %>% stringr::str_replace("^", "Root;")

# read in sequences passing max thresholds
seqs_max <- Biostrings::readDNAStringSet(seqs_max_fasta)

# consensus parameters
con_min_n <- as.numeric(con_min_n)
con_min_prop <- as.numeric(con_min_prop)

### run code

allowed_ranks <- c("kingdom","phylum", "class", "order", "family", "genus", "species")
root_ranks <- c("root", allowed_ranks)

# get taxonomy of consistent genera
cg_tax <- 
	cg_all %>%
	# get family, order, class, phylum and kingdom
	dplyr::mutate(
		genus = paste0("Root;",genus),
		family = stringr::str_extract(genus, ".+(?=;)"),
        order = stringr::str_extract(family, ".+(?=;)"),
        class = stringr::str_extract(order, ".+(?=;)"),
        phylum = stringr::str_extract(class, ".+(?=;)"),
        kingdom = stringr::str_extract(phylum, ".+(?=;)"),
		root = "Root"
	)

# get consistent genera that are "large" (>= 2 * con_min_n representative sequences) 
lcg_tax <- 
	cg_tax %>%
	dplyr::filter(rep_n >= 2*con_min_n) 

### loop through each consistent genus, finding up to 3 random genera per LSR

root.to.genus <- root_ranks[-8]

comparator_list <- 
	lapply(
		seq_along(cg_vec),
		function(y){
			x <- cg_vec[y]
			# get central sequence of focal genus
			x_central <- cg_tax[cg_tax$genus == x, "name"] %>% as.character()
			# get vector of focal genus taxa at each rank, with "Unclassified" as "UNCLASSIFIED"
			x_tax <- 
				stringr::str_split(x, ";", n = 7, simplify = F) %>% 
			    unlist() %>% 
			    stringr::str_replace_all(., "^Unclassified$", "UNCLASSIFIED")
			# which ranks are classified?
			x_classified <- 
			    x_tax %>% 
				stringr::str_detect(., "UNCLASSIFIED|Unclassified", negate = T)
			# get lowest classified rank for focal sequence (ie. max index)
			x_lcr_i <- which(root.to.genus == root.to.genus[max(match(root.to.genus[x_classified], root.to.genus))] )
			# classify all LCG by LSR with focal genus
			selections <- 
				lcg_tax %>%
				# remove focal genus
				dplyr::filter(!genus == x) %>%
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
				# select up to 3 genera per rank
				dplyr::slice_sample(n = 3)
			# combine selected sequences with focal central sequence
			out <- c(x_central,selections$name)
			if (y %% 100 == 0) message(stringr::str_glue("Selected comparators for {y} of {length(cg_vec)} genera"))
			# if no comparators were found, return nothing
			if (length(out) < 2){
				return(NULL)
			} else {
				return(out)
			}
		}
	)

# remove null elements
comparator_list[sapply(comparator_list, is.null)] <- NULL

message("Converting selection vectors to .cfa format")

# convert each comparator vector into a DSS object then convert to .cfa format
if (length(comparator_list) > 0){
	cfa_out <- 
		lapply(
			comparator_list,
			function(x){
				# use names for subsetting to get sequences in order of the vector (ie. focal first, then comparators)
				x_seqs_char <- seqs_max[x] %>% as.character()
				x_seqs_cfa <- 
					as.vector(rbind(names(x_seqs_char) %>% stringr::str_replace("^",">"),unname(x_seqs_char))) %>% 
					base::paste(collapse = ">>>")
				return(x_seqs_cfa)
			}
		) %>%
		unlist()
} else {
	cfa_out <- NULL
}

# writes empty file if needed
readr::write_lines(cfa_out, "min_comparators.cfa")