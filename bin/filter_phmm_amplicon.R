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
    "primer_info_file",
    "hmmer_output",
    "trim_to_amplicon",
    "amplicon_min_length",
    "amplicon_min_cov",
    "remove_primers"
)
lapply(nf_vars, nf_var_check)

### process variables 
# read in fasta file as DNAbin
seqs <- ape::read.FASTA(fasta_file)

# read in hmmer_output
domtblout <- readr::read_lines(hmmer_output, lazy=FALSE, progress=FALSE) 

# read in frame info file
primer_info <- readr::read_csv(primer_info_file)

# filtering parameters
amplicon_min_length <- amplicon_min_length %>% as.integer()
amplicon_min_cov <- amplicon_min_cov %>% as.numeric()
if(trim_to_amplicon == "true"){ trim_to_amplicon <- TRUE } else if(trim_to_amplicon == "false") { trim_to_amplicon <- FALSE } else {stop(paste0("'--trim_to_amplicon' must be 'true' or 'false'"))}
if(remove_primers == "true"){ remove_primers <- TRUE } else if(remove_primers == "false") { remove_primers <- FALSE } else {stop(paste0("'--remove_primers' must be 'true' or 'false'"))}

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
# code inspired by rhmmer: https://github.com/arendsee/rhmmer
domtblout_parsed <- 
    domtblout %>%
    # remove all rows that begin with # as they aren't part of the table
    .[stringr::str_detect(., "^#", negate = TRUE)] %>% 
    # replace runs of spaces with '\t'
    stringr::str_replace_all(., "\\s+", "\t") %>% 
    # collapse vector to a single element
    paste0(collapse = "\n") %>%
    # add new line character to end of element
    stringr::str_replace(., "$", "\n") %>%
    # read in element as if it were a tsv file
    readr::read_tsv(
        col_names = names(col_types$cols),
        na = '-',
        col_types = col_types,
        lazy = FALSE, 
        progress = FALSE
    ) %>% 
    # combine additional end columns with "description" (if multiple instances of space characters in species name)
    tidyr::unite(description, c(description, dplyr::starts_with("X")), sep = " ") %>%
    dplyr::mutate(description = dplyr::if_else(description == "NA", NA_character_, description)) %>%
    # replace '\t' in description column with spaces to regenerate sequence names
    dplyr::mutate(
        description = stringr::str_replace_all(description, "\\\t", " "),
        target_name = stringr::str_replace_all(target_name, "\\!\\?\\!\\?", " ") # replace each "!?!?" with a space character as well 
    ) 

## calculate the expected size of the specified amplicon based on PHMM size, start/end offsets, primer padding and primer removal
# NOTE: This might be smaller than the actual amplicon from most sequences, due to indels relative to the PHMM
hmm_length <- domtblout_parsed$query_len %>% unique() %>% {.*3}
if (remove_primers == TRUE){
    amplicon_length <- 
        hmm_length - primer_info$offset_start + primer_info$offset_end - primer_info$plen_start - primer_info$plen_end
} else {
    amplicon_length <- 
        hmm_length - primer_info$offset_start + primer_info$offset_end + primer_info$pad_start + primer_info$pad_end
}

# get sequence lengths in bases
seq_lengths <- 
    lengths(seqs) %>% 
    tibble::enframe(name = "target_name", value = "bases")

# convert to single sequence hits
hits <- 
    domtblout_parsed %>%
    # undo HMMER splitting sequence names into two parts at the first space character
    dplyr::mutate(
        target_name = dplyr::if_else(
            !is.na(description),
            paste0(target_name, " ", description), 
            target_name
        )
    ) %>%
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
    ) %>%
    # join sequence nucleotide lengths
    dplyr::left_join(., seq_lengths, by = "target_name") %>%
    # convert envelope coordinates to nucleotide coordinates
    dplyr::mutate(
        ## force alignment to the ends of the PHMM, if possible
        # subtract start gap from alignment start position, not going lower than start of translated sequence
        force_from = base::pmax(ali_from - hmm_start_gap, 1),
        # add end gap to alignment end position, not going higher than end of translated sequence
        force_to = base::pmin(ali_to + hmm_end_gap, target_len),
        # recalculate the start gap by seeing how much the start of the alignment moved and subtracting from original start gap
        hmm_start_gap_new = hmm_start_gap - (ali_from - force_from),
        # recalculate the end gap by seeing how much the end of the alignment moved and subtracting from original end gap
        hmm_end_gap_new = hmm_end_gap - (force_to - ali_to)
    ) %>%
    # calculate nucleotide positions of the forced alignment, taking primers into consideration
    { if(remove_primers == TRUE){
        dplyr::mutate(
            .data = .,
            nuc_from = 
                # don't go lower than 1
                base::pmax(
                    1,
                    # nucleotide position of the start of the forced aa alignment
                    (force_from * 3) - 2 +
                        # add (to position) remainder of primer length after subtracting hmm start gap plus start offset 
                        base::pmax(0, primer_info$plen_start + primer_info$offset_start - (hmm_start_gap_new * 3))
                ),
            nuc_to =
                # don't go higher than nucleotide length 
                base::pmin(
                    bases,
                    # nucleotide position of the start of the forced aa alignment
                    (force_to * 3)  -
                        # subtract (from position) remainder of primer length after subtracting hmm end gap plus end offset
                        base::pmax(0, primer_info$plen_end - primer_info$offset_end - (hmm_end_gap_new * 3))
                ),
            nuc_len = (nuc_to - nuc_from) + 1
        )
    } else {
        dplyr::mutate(
            .data = .,
            nuc_from = 
                # don't go lower than 1
                base::pmax(
                    1,
                    # nucleotide position of the start of the forced aa alignment
                    ((force_from * 3) - 2) + 
                        # only add start offset and padding if larger than hmm start gap in bp
                        # note: start offset is negative
                        base::pmin(0, primer_info$offset_start + (hmm_start_gap_new * 3) - primer_info$pad_start)
                ),
            nuc_to =
                # don't go higher than nucleotide length 
                base::pmin(
                    bases,
                    # nucleotide position of the end of the forced aa alignment
                    (force_to * 3) +
                        # only add end offset and padding if larger than hmm end gap in bp
                        base::pmax(0, primer_info$offset_end - (hmm_end_gap_new * 3) + primer_info$pad_end)
                ),
            nuc_len = (nuc_to - nuc_from) + 1
        )
    } } %>%
    dplyr::mutate(
      # coverage of record amplicon vs idealised amplicon
      amp_cov = nuc_len / amplicon_length,
      # coverage of record amplicon hit vs PHMM 
      phmm_cov = (hmm_to - hmm_from + 1) / query_len
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
       nuc_len >= amplicon_min_length &
       phmm_cov >= amplicon_min_cov
    ) 

# subset to retained sequences
if (any(names(seqs) %in% hits_retained$target_name)){
    seqs_combined <- seqs[names(seqs) %in% hits_retained$target_name]
} else {
    seqs_combined <- list() %>% as.DNAbin()
}

hit_locations <- 
        hits_retained %>%
        # reorder tibble to match order of combined DNAbin
        dplyr::arrange(
            factor(target_name, levels = names(seqs_combined))
        ) 

# trim nucleotide seqs based on location of hits for each sequence
seqs_subset <-
    seqs_combined %>%
    DNAbin2DNAstringset(.) %>%
    XVector::subseq(., start = hit_locations$nuc_from, end = hit_locations$nuc_to) %>%
    ape::as.DNAbin()

# check calculated subset lengths are the same as the ones done
if (!all(unname(lengths(seqs_subset)) == hit_locations$nuc_len)){
    stop("Actual subset sequence lengths are not the same as calculated subset sequence lengths")
}

# seqs without a significant hit 
seqs_nohit <- seqs[!names(seqs) %in% hits$target_name]

# save data for hits that were removed after filtering
seqs_removed_tibble <- 
    hits %>%
    # keep only sequences with failed hits
    dplyr::filter(
       nuc_len < amplicon_min_length |
       phmm_cov < amplicon_min_cov
    ) %>%
    # which filters were failed
    dplyr::mutate(
        fail_min_length = nuc_len < amplicon_min_length,
        fail_min_cov = phmm_cov < amplicon_min_cov
    ) %>%
    dplyr::mutate(removal_type = "excluded_hit", .after = phmm_cov) %>%
    # add names of sequences without a hit
    tibble::add_row(
        target_name = names(seqs_nohit)[!names(seqs_nohit) %in% .$target_name], 
        removal_type = "no_hit",
    ) %>%
    dplyr::select(-old_name) # remove old_name as redundant with target_name

readr::write_csv(seqs_removed_tibble, file = "removed_amplicon.csv")

### outputs

# removed seqs
seqs_removed <- seqs[names(seqs) %in% seqs_removed_tibble$target_name]

# check all sequences are either retained or removed
if ( length(seqs_subset) + nrow(seqs_removed_tibble) != length(seqs) ){
    stop("ERROR: Sum of retained and removed sequences does not equal the number of input sequences")
}

# write fasta of subsetted sequences
if ( !is.null(seqs_subset) && length(seqs_subset) > 0 ){
    write_fasta(
        seqs_subset, 
        file = paste0("retained_amplicon.fasta"), 
        compress = FALSE
        )
} else {
    file.create(paste0("retained_amplicon.fasta"))
}

# write fasta of excluded (filtered-out) sequences
if ( !is.null(seqs_removed) && length(seqs_removed) > 0 ){
    write_fasta(
        seqs_removed, 
        file = paste0("removed_amplicon.fasta"), 
        compress = FALSE
        )
} else {
    file.create(paste0("removed_amplicon.fasta"))
}