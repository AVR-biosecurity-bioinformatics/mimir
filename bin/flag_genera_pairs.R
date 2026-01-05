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
    "alignment_file",
    "thresholds_file",
    "seqs_file",
    "counts_file"
)
lapply(nf_vars, nf_var_check)

### process variables 

# read in list of sequences
alignment_single <- readr::read_lines(alignment_file)

# process single-line alignment format into DSS format
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

thresholds <- readr::read_csv(thresholds_file, show_col_types = F)

seqs <- Biostrings::readDNAStringSet(seqs_file)

counts <- readr::read_tsv(counts_file, col_names = c("name","n"), show_col_types = FALSE)

allowed_ranks <- c("kingdom","phylum", "class", "order", "family", "genus", "species")
root_ranks <- c("root", allowed_ranks)

### run code

# for each DSS object in seqs_list, produce a distance matrix, then pull the distance from the first sequence (which is the query)
flags <-
	lapply(
		seq_along(seqs_list),
		function(x){
			x_seqs <- seqs_list[[x]]
			if (length(x_seqs) > 1){
				# get distance matrix
				x_d <- DECIPHER::DistanceMatrix(x_seqs, verbose = F)
				# parse taxonomy of query
				query_name <- names(x_seqs)[1]
				query_tax <- 
				    stringr::str_split(query_name, ";", n = 8, simplify = F) %>% 
				    unlist %>% 
				    stringr::str_replace_all(., "^Unclassified$", "UNCLASSIFIED") %>% 
				    .[2:8] %>%
				    c("Root", .)
				query_k <- base::paste(query_tax[1:2], collapse = ";")
				query_p <- base::paste(query_tax[1:3], collapse = ";")
				query_c <- base::paste(query_tax[1:4], collapse = ";")
				query_o <- base::paste(query_tax[1:5], collapse = ";")
				query_f <- base::paste(query_tax[1:6], collapse = ";")
				query_g <- base::paste(query_tax[1:7], collapse = ";")
				query_s <- base::paste(query_tax[1:8], collapse = ";")
				query_classified <- stringr::str_detect(query_tax, "UNCLASSIFIED|Unclassified(\\d+)?", negate = T)
				# get lowest classified rank for focal sequence (ie. max index)
				query_lcr <- root_ranks[max(match(root_ranks[query_classified], root_ranks))] 
				query_lcr_i <- which(root_ranks == query_lcr)
				# find pairs of sequences that violate the threshold for their LSR
				x_flags <- 
				    x_d[1,] %>%
				    tibble::enframe(name = "target", value = "fid") %>%
				    dplyr::mutate(
				        fid = 1 - fid,
				        query = query_name
				    ) %>%
				    # keep only top 100 hits, not including self
				    dplyr::filter(target != query) %>%
				    # extract taxonomy of target sequences
				    dplyr::mutate(
				        species = paste("Root", stringr::str_extract(target, "(?<=;).+$"), sep = ";"),
				        genus = stringr::str_extract(species, ".+(?=;)"),
				        family = stringr::str_extract(genus, ".+(?=;)"),
				        order = stringr::str_extract(family, ".+(?=;)"),
				        class = stringr::str_extract(order, ".+(?=;)"),
				        phylum = stringr::str_extract(class, ".+(?=;)"),
				        kingdom = stringr::str_extract(phylum, ".+(?=;)"),
				        lsr = dplyr::case_when(
				            query_s == species ~ "species",
				            query_g == genus ~ "genus",
				            query_f == family ~ "family",
				            query_o == order ~ "order",
				            query_c == class ~ "class",
				            query_p == phylum ~ "phylum",
				            query_k == kingdom ~ "kingdom",
				            .default = "root"
				        ),
				        # lowest classified rank for the target sequence
				        lcr = dplyr::case_when(
				            stringr::str_detect(species, "Unclassified", negate = T) ~ "species",
				            stringr::str_detect(genus, "Unclassified", negate = T) ~ "genus",
				            stringr::str_detect(family, "Unclassified", negate = T) ~ "family",
				            stringr::str_detect(order, "Unclassified", negate = T) ~ "order",
				            stringr::str_detect(class, "Unclassified", negate = T) ~ "class",
				            stringr::str_detect(phylum, "Unclassified", negate = T) ~ "phylum",
				            stringr::str_detect(kingdom, "Unclassified", negate = T) ~ "kingdom",
				            .default = "root"
				        )
				    ) %>%
				    dplyr::rowwise() %>%
				    dplyr::mutate(
				        # get the lowest rank that is classified for both sequences
				        lcr = root_ranks[min(query_lcr_i, which(root_ranks == lcr))]
				    ) %>%
				    dplyr::ungroup() %>%
				    dplyr::select(query, target, lsr, lcr, fid) %>%
				    dplyr::mutate(
				        # is fid above max for the lsr rank?
				        threshold = dplyr::case_when(
				            lsr == "family" & lcr %in% root_ranks[7:8] ~ as.numeric(thresholds[thresholds$lsr == "family", "max"]),
				            lsr == "order" & lcr %in% root_ranks[6:8] ~ as.numeric(thresholds[thresholds$lsr == "order", "max"]),
				            lsr == "class" & lcr %in% root_ranks[5:8] ~ as.numeric(thresholds[thresholds$lsr == "class", "max"]),
				            lsr == "phylum" & lcr %in% root_ranks[4:8] ~ as.numeric(thresholds[thresholds$lsr == "phylum", "max"]),
				            lsr == "kingdom" & lcr %in% root_ranks[3:8] ~ as.numeric(thresholds[thresholds$lsr == "kingdom", "max"]),
				            lsr == "root" & lcr %in% root_ranks[2:8] ~ NA,
				            .default = NA
				        ),
				        flag = fid > threshold
				    ) %>% 
				    dplyr::filter(flag == TRUE)
			} else {
				x_flags <-
					tibble::tibble(
						query = character(),
						target = character(),
						lsr = character(),
						lcr = character(),
						fid = numeric(),
						threshold = numeric(),
						flag = logical()
					)
			}
			if (x %% 10 == 0){
				message(stringr::str_glue("Finished {x} of {length(seqs_list)}"))
			}
			return(x_flags)
		}
	) %>%
	dplyr::bind_rows()

# convert flags to genera pairs
genera_pairs <- 
	flags %>%
	# for each query and LSR, get the target with the highest identity
	dplyr::group_by(query, lsr) %>%
	dplyr::arrange(desc(fid), .by_group = T) %>%
	dplyr::slice(1) %>%
	dplyr::ungroup() %>%
	dplyr::mutate(
		g_q = stringr::str_extract(query, "(?<=;).+$") %>% stringr::str_extract(., ".+(?=;)"), 
		g_t = stringr::str_extract(target, "(?<=;).+$") %>% stringr::str_extract(., ".+(?=;)")
	) %>%
	dplyr::select(query = g_q, target = g_t) %>%
	dplyr::distinct()

# write out
readr::write_csv(genera_pairs, "genera_pairs.csv")
