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
    "bold_seqs_file",
    "genbank_seqs_file",
    "rankedlineage_noname",
    "gb_dups_file"
    )
lapply(nf_vars, nf_var_check)

### process variables 

bold_seqs <- Biostrings::readDNAStringSet(bold_seqs_file)

genbank_seqs <- Biostrings::readDNAStringSet(genbank_seqs_file)

gb_dups <- readr::read_csv(gb_dups_file)

rankedlineage_noname <- readRDS(rankedlineage_noname)

if (!prefer_source %in% c("bold","genbank")) stop("'--prefer_source' must be 'bold' or 'genbank'")
if (!duplicate_taxonomy %in% c("keep_source","lca")) stop("'--duplicate_taxonomy' must be 'keep_source' or 'lca'")
if (!duplicate_sequence %in% c("ignore_sequence","longest","unambiguous_n","unambiguous_prop")) stop("'--duplicate_sequence' must be 'ignore_sequence', 'longest', 'unambiguous_n' or 'unambiguous_prop'")

### custom functions
# function for getting lowest common ancestor of two lineages
lin_lca <- function(x, y){
	x_split <- stringr::str_split_1(x, ";") 
	y_split <- stringr::str_split_1(y, ";") %>% stringr::str_replace(., "Unclassified", "UNCLASSIFIED")
	out <- replace(x_split, which(x_split != y_split), values = "Unclassified") %>% stringr::str_flatten(., collapse = ";")
	return(out)
}

# function for selecting the source from sequence
select_sequence_source <- function(sb, sg, type, prefer = prefer_source){
	if (!prefer %in% c("bold","genbank")) stop("Incorrect prefer argument")
	if (type == "ignore_sequence"){
		return(prefer)
	} else if (type == "longest"){
		sb_len <- stringr::str_length(sb)
		sg_len <- stringr::str_length(sg)
		
		if (sb_len > sg_len)			return("bold")
		else if (sg_len > sb_len)		return("genbank")
		else if (sb_len == sg_len)		return(prefer)
	
	} else if (type == "unambiguous_n"){
		n_sb <- stringr::str_count(sb, "[ACTG]")
		n_sg <- stringr::str_count(sg, "[ACTG]")
		
		if (n_sb > n_sg)				return("bold")
		else if (n_sg > n_sb)			return("genbank")
		else if (n_sg == n_sb)			return(prefer)
		
	} else if (type == "unambiguous_prop"){
		prop_sb <- stringr::str_count(sb, "[ACTG]") / stringr::str_length(sb)
		prop_sg <- stringr::str_count(sg, "[ACTG]") / stringr::str_length(sg)
	
		if (prop_sb > prop_sg)			return("bold")
		else if (prop_sg > prop_sb) 	return("genbank")
		else if (prop_sg == prop_sb)	return(prefer)
		
	} else stop("Incorrect type argument") 
}

### run code

# get a tibble of NCBI taxids for all lineages
taxid_lin <- 
	rankedlineage_noname %>%
	dplyr::mutate(
        dplyr::across(species:kingdom, .fns = ~replace(., is.na(.), "Unclassified")),
        lin_new = paste0(kingdom,";",phylum,";",class,";",order,";",family,";",genus,";",species)
    ) %>%
	dplyr::select(taxid_new = tax_id, lin_new) %>%
	dplyr::arrange(taxid_new) %>% 
	dplyr::distinct(lin_new, .keep_all = T)

## create tibbles of bold and genbank sequence lineages with sequences
bold_tibble <- 
	bold_seqs %>%
	as.character() %>%
	tibble::enframe(name = "name_bold", value = "sequence_bold") %>%
	dplyr::mutate(
		seqid_bold = stringr::str_extract(name_bold, "^.+(?=\\|)"),
		taxid_bold = stringr::str_extract(name_bold, "(?<=\\|).+?(?=;)"),
		lin_bold = stringr::str_extract(name_bold, "(?<=;).+$")
	) %>%
	dplyr::distinct() # remove duplicates if there are any
	
genbank_tibble <- 
	genbank_seqs %>% 
	as.character() %>%
	tibble::enframe(name = "name_genbank", value = "sequence_genbank") %>%
	dplyr::mutate(
		seqid_genbank = stringr::str_extract(name_genbank, "^.+(?=\\|)"),
		taxid_genbank = stringr::str_extract(name_genbank, "(?<=\\|).+?(?=;)"),
		lin_genbank = stringr::str_extract(name_genbank, "(?<=;).+$"),
		acc_genbank = stringr::str_extract(name_genbank, "^.+(?=\\.)")
	) %>%
	dplyr::distinct() # remove duplicates if there are any

if (genbank_tibble$acc_genbank %>% duplicated() %>% any()){
	stop("Some GenBank accessions are duplicated")
}

## match up bold sequences that have GB duplicates with their duplicate sequences
dups_joined <- 
	gb_dups %>%
	dplyr::rename(name_bold = name, acc_genbank = genbank_acc) %>%
	dplyr::left_join(., genbank_tibble, by = "acc_genbank") %>%
	dplyr::select(-acc_genbank) %>%
	dplyr::left_join(., bold_tibble, by = "name_bold") %>%
	dplyr::mutate(
		lin_match = lin_bold == lin_genbank,
		seq_match = sequence_bold == sequence_genbank
	)

# work out which duplicate record to keep, and how taxonomy is kept or changed
# evaluate: 
# 1. duplicate_sequence to select source
# 2. in case of 'ignore_sequence' or tie at 1., use prefer_source to select source
# 3. change taxonomy if duplicate_taxonomy == 'lca', otherwise use source taxonomy
dups_processed <- 
	dups_joined %>%
	# group rowwise as custom functions aren't vectorised
	dplyr::rowwise() %>%
	dplyr::mutate(
		source_new = select_sequence_source(
			sequence_bold, 
			sequence_genbank, 
			type = duplicate_sequence, 
			prefer = prefer_source
		),
		lin_new = 
			if(duplicate_taxonomy == "lca")     lin_lca(lin_bold, lin_genbank)
			else if(source_new == "bold")       lin_bold
			else if(source_new == "genbank")    lin_genbank
			else                                NA
	) %>%
	dplyr::ungroup() %>%
	dplyr::mutate(
		seqid_new = dplyr::case_match(
			source_new, 
			"bold" ~ seqid_bold,
			"genbank" ~ seqid_genbank,
			.default = NA
		),
		sequence_new = dplyr::case_match(
			source_new, 
			"bold" ~ sequence_bold,
			"genbank" ~ sequence_genbank,
			.default = NA
		)
	) %>%
	dplyr::select(source_new, name_bold, name_genbank, taxid_bold, taxid_genbank, seqid_new, lin_new, sequence_new) %>%
	# match NCBI taxids to taxid_lin
	dplyr::left_join(., taxid_lin, by = "lin_new") %>%
	# if an NCBI taxid matches, use that, otherwise keep old taxid
	dplyr::mutate(
		taxid_new = dplyr::case_when(
			!is.na(taxid_new) ~ paste0("NCBI:",taxid_new),
			source_new == "bold" ~ taxid_bold,
			source_new == "genbank" ~ taxid_genbank,
			.default = NA
		),
		# construct new record name
		name_new = paste0(seqid_new,"|",taxid_new,";",lin_new)
	) %>%
	dplyr::select(source_new, seqid_new, taxid_new, lin_new, sequence_new, name_bold, name_genbank, name_new)

# if any new fields are NA, throw error
dups_na_check <- 
	dups_processed %>%
	dplyr::filter(is.na(seqid_new) | is.na(sequence_new) | is.na(taxid_new) | is.na(lin_new) )

if (nrow(dups_na_check) > 0) stop("Duplicate processing produced one or more NA values for new fields")
	
# split sequences into BOLD source and GenBank source as DSS objects
processed_bold <- 
	dups_processed %>% 
	dplyr::filter(source_new == "bold") %>%
	dplyr::select(name_new, sequence_new) %>%
	tibble::deframe() %>% 
	Biostrings::DNAStringSet(.)

processed_genbank <- 
	dups_processed %>% 
	dplyr::filter(source_new == "genbank") %>%
	dplyr::select(name_new, sequence_new) %>%
	tibble::deframe() %>% 
	Biostrings::DNAStringSet(.)

# add processed duplicates (per source) to input sequences without the duplicates
bold_seqs_out <- 
    c(
        bold_seqs[!names(bold_seqs) %in% dups_processed$name_bold], 
        processed_bold
    )

genbank_seqs_out <- 
    c(
        genbank_seqs[!names(genbank_seqs) %in% dups_processed$name_genbank], 
        processed_genbank
    )

# write sequences to file
Biostrings::writeXStringSet(bold_seqs_out, "bold_out.fasta", width = 9999)
Biostrings::writeXStringSet(genbank_seqs_out, "genbank_out.fasta", width = 9999)

# export processing tibble for debugging
readr::write_csv(dups_processed, "bold_gb_duplicates.csv")

gc()