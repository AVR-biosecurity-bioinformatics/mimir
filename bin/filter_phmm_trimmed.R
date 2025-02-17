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
    "fasta_file",
    "hmmer_output",
    "hmm_max_evalue", 
    "hmm_min_score", 
    "hmm_max_hits",   
    "hmm_min_acc",   
    "hmm_max_gap",
    "min_length_trimmed"
)
lapply(nf_vars, nf_var_check)

### process variables 
# read in fasta file as DNAbin
seqs <- ape::read.FASTA(fasta_file)

# read in hmmer_output
domtblout <- readr::read_lines(hmmer_output, lazy=FALSE, progress=FALSE) 

# filtering parameters
evalue_threshold <- hmm_max_evalue %>% as.numeric()
score_threshold <- hmm_min_score %>% as.numeric()
domain_threshold <- hmm_max_hits %>% as.integer() # max allowed domains (1 is safe)
acc_threshold <- hmm_min_acc %>% as.numeric()
terminalgap_threshold <- hmm_max_gap %>% as.integer() # terminal gap is the distance between the beginning or end of the envelope/HMM and the min/max of the target or HMM
min_length_trimmed <- min_length_trimmed %>% as.integer()
## NOTE: terminal gaps above the threshold on the same side (start or end) for both the envelope and HMM is expected in the case of frameshifts
## More extensive frameshift detection could be done with alignments to other sequences and removing those with gaps not in multiples of three

# alter filtering parameters for trimmed
evalue_mantissa <- evalue_threshold %>% as.character() %>% stringr::str_extract(., "^\\d+") %>% as.double
evalue_power <- evalue_threshold %>% as.character() %>% stringr::str_extract(., "\\d+$") %>% as.integer 
evalue_threshold_trimmed <- paste0(evalue_mantissa,"e-",(evalue_power/2)) %>% as.numeric()
if (evalue_threshold_trimmed >1){stop("New evalue is greater than 1")}

score_threshold_trimmed <- floor(score_threshold / 5)
terminalgap_threshold_trimmed <- ceiling(terminalgap_threshold / 2)

### run code

## import and parse HMMER output
# create column specifications
col_types <- 
    readr::cols(
        target_name         = readr::col_character(), # name of target seq
        target_accession    = readr::col_character(), # accession of target (unused)
        target_len          = readr::col_double(), # length of target (in residues)
        query_name          = readr::col_character(), # name of query profile
        query_accession     = readr::col_character(), # accession of query profile
        query_len           = readr::col_double(), # length of query profile (in residues)
        evalue_full         = readr::col_double(), # e-value of full target sequence (across all domains)
        score_full          = readr::col_double(), # score of full target sequence (across all domains)
        bias_full           = readr::col_double(), # bias of full target sequence (across all domains)
        domain_n            = readr::col_integer(), # this domain's number
        domain_total        = readr::col_integer(), # total number of domains in the target sequence
        evalue_c            = readr::col_double(), # conditional e-value of domain (lenient)
        evalue_i            = readr::col_double(), # independent e-value of domain (stringent)
        score               = readr::col_double(), # score for domain
        bias                = readr::col_double(), # bias for domain
        hmm_from            = readr::col_double(), # position in hmm profile the hit starts
        hmm_to              = readr::col_double(), # ... ends
        ali_from            = readr::col_double(), # position in the target the hit starts
        ali_to              = readr::col_double(), # ... ends
        env_from            = readr::col_double(), # position in target when envelope (fuzzy match) starts
        env_to              = readr::col_double(), # ... ends
        acc                 = readr::col_double(), # mean posterior probability of aligned residues in the MEA alignment; 1.00 = completely reliable
        description         = readr::col_character() # free text description
    )

N <- length(col_types$cols)

# parse space-delimited hmmer output
# code from rhmmer: https://github.com/arendsee/rhmmer
domtblout_parsed <- 
    domtblout %>%
    sub(
        pattern = sprintf("(%s).*", paste0(rep('\\S+', N), collapse=" +")),
        replacement = '\\1',
        x = .,
        perl = TRUE
    ) %>%
    gsub(pattern = "  *", replacement = "\t") %>%
    paste0(collapse = "\n") %>%
    readr::read_tsv(
        col_names = names(col_types$cols),
        comment = '#',
        na = '-',
        col_types = col_types,
        lazy = FALSE,
        progress = FALSE
    )

# convert to single sequence hits
hits <- 
    domtblout_parsed %>%
    # rename 'target_name' to retain original name for later
    dplyr::rename( old_name = target_name ) %>%
    # extract stop codons from target_name
    tidyr::separate(
        col = old_name, 
        into = c("target_name","stop_codons"), 
        sep = "_sc=", 
        remove = FALSE
    ) %>% 
    # split target name into original sequence name and the frame number 
    tidyr::separate(
        col = target_name, 
        into = c("target_name","frame"), 
        sep = "_frame=", 
        remove = TRUE
    ) %>% 
    # reformat new columns
    dplyr::mutate(
        frame = as.integer(frame),
        stop_codons = dplyr::na_if(stop_codons, "") %>% as.character() # empty value to NA
        ) %>%
    # group by sequence, arranging from best score to worst score (if there are multiple hits per frame), keeping best hit
    dplyr::group_by(target_name) %>%
    dplyr::arrange(desc(score)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    # calculate envelope and hmm gaps
    dplyr::mutate(
        env_start_gap = env_from - 1,
        env_end_gap = target_len - env_to,
        hmm_start_gap = hmm_from - 1,
        hmm_end_gap = query_len - hmm_to
    ) 

# check hits are all in sequence file
if(!all(hits$target_name %in% names(seqs))){
    stop("One or more sequence names in HMMER output table do not match names in input .fasta file")
}


# filter hits based on various thresholds
hits_retained <- 
    hits %>%
    # remove sequences that don't meet all criteria
    dplyr::filter(
        evalue_i <= evalue_threshold_trimmed & 
        score >= score_threshold_trimmed & 
        # domain_total == domain_threshold & 
        acc >= acc_threshold &
        !( hmm_start_gap > terminalgap_threshold_trimmed & env_start_gap > terminalgap_threshold_trimmed ) & # both need to be fulfilled 
        !( hmm_end_gap > terminalgap_threshold_trimmed & env_end_gap > terminalgap_threshold_trimmed ) 
    )

# fwd seqs
if (any(names(seqs) %in% hits_retained$target_name)){
    seqs_combined <- seqs[names(seqs) %in% hits_retained$target_name]
} else {
    seqs_combined <- list() %>% as.DNAbin()
}

# get sequence lengths in bases
seq_lengths <- 
    lengths(seqs_combined) %>% 
    tibble::enframe(name = "target_name", value = "bases")

# calculate positions of hits on nucleotide sequence
hit_locations <- 
    hits_retained %>%
    # reorder tibble to match order of combined DNAbin
    dplyr::arrange(
        factor(target_name, levels = names(seqs_combined))
    ) %>%
    # join sequence nucleotide lengths
    dplyr::left_join(., seq_lengths, by = "target_name") %>%
    # convert envelope coordinates to nucleotide coordinates
    dplyr::mutate(
        nuc_from =  (env_from * 3) - (3 - abs(frame)),  
        nuc_to =    (env_to * 3) + (abs(frame) - 1),
        nuc_len =   nuc_to - nuc_from + 1,
        coverage =  nuc_len / bases,
        pad_from =  base::pmax(1,nuc_from-2), # pad hit location by up to 2 bases 5'
        pad_to =    base::pmin(bases,nuc_to+2), # pad hit location by up to 2 bases 3'
        pad_len =   pad_to - pad_from + 1
    ) 
  

# subset nucleotide seqs based on location of hits for each sequence
seqs_subset <-
    seqs_combined %>%
    DNAbin2DNAstringset(.) %>%
    XVector::subseq(., start = hit_locations$pad_from, end = hit_locations$pad_to) %>%
    ape::as.DNAbin()

# check calculated subset lengths are the same as the ones done
if (!all(unname(lengths(seqs_subset)) == hit_locations$pad_len)){
    stop("Actual subset sequence lengths are not the same as calculated subset sequence lengths")
}

# remove sequences shorter than min_length_trimmed
seqs_long <- seqs_subset[seqs_subset %>% lengths() %>% unname() >= min_length_trimmed]

# seqs without a significant hit (either no hit or an excluded hit)
seqs_nohit <- seqs[!names(seqs) %in% hits_retained$target_name]

# save data for hits that were removed after filtering
seqs_removed_tibble <- 
    hits %>%
    # remove sequences that don't meet all criteria
    dplyr::filter(
        evalue_i > evalue_threshold_trimmed |
        score < score_threshold_trimmed | 
        # domain_total > domain_threshold |
        acc < acc_threshold |
        ( hmm_start_gap > terminalgap_threshold_trimmed & env_start_gap > terminalgap_threshold_trimmed ) | # both need to be fulfilled 
        ( hmm_end_gap > terminalgap_threshold_trimmed & env_end_gap > terminalgap_threshold_trimmed ) 
    ) %>%
    # which filters were failed
    dplyr::mutate(
        fail_max_evalue =       evalue_i > evalue_threshold_trimmed,
        fail_min_score =        score < score_threshold_trimmed,
        # fail_max_hits =         domain_total > domain_threshold ,
        fail_min_acc =          acc < acc_threshold ,
        fail_max_gap_start =    ( hmm_start_gap > terminalgap_threshold_trimmed & env_start_gap > terminalgap_threshold_trimmed ),
        fail_max_gap_end =      ( hmm_end_gap > terminalgap_threshold_trimmed & env_end_gap > terminalgap_threshold_trimmed )
    ) %>%
    dplyr::mutate(removal_type = "excluded_hit", .after = hmm_end_gap) %>%
    # add names of sequences without a hit
    tibble::add_row(
        target_name = names(seqs_nohit)[!names(seqs_nohit) %in% .$target_name], 
        removal_type = "no_hit",
    ) %>%
    # add sequences removed for being too short
    tibble::add_row(
        target_name = names(seqs_subset[seqs_subset %>% lengths() %>% unname() < min_length_trimmed]),
        removal_type = "length",
        fail_max_evalue =       FALSE,
        fail_min_score =        FALSE,
        # fail_max_hits =         FALSE,
        fail_min_acc =          FALSE,
        fail_max_gap_start =    FALSE,
        fail_max_gap_end =      FALSE,
    ) %>%
    dplyr::select(-old_name) # remove old_name as redundant with target_name

readr::write_csv(seqs_removed_tibble, file = "removed_trimmed.csv")

### outputs

# removed seqs
seqs_removed <- seqs[names(seqs) %in% seqs_removed_tibble$target_name]

# check all sequences are either retained or removed
if ( length(seqs_long) + nrow(seqs_removed_tibble) != length(seqs) ){
    stop("ERROR: Sum of retained and removed sequences does not equal the number of input sequences")
}

# write fasta of subsetted sequences
if ( !is.null(seqs_long) && length(seqs_long) > 0 ){
    write_fasta(
        seqs_long, 
        file = paste0("retained_trimmed.fasta"), 
        compress = FALSE
        )
} else {
    file.create(paste0("retained_trimmed.fasta"))
}

# write fasta of excluded (filtered-out) sequences
if ( !is.null(seqs_removed) && length(seqs_removed) > 0 ){
    write_fasta(
        seqs_removed, 
        file = paste0("removed_trimmed.fasta"), 
        compress = FALSE
        )
} else {
    file.create(paste0("removed_trimmed.fasta"))
}
