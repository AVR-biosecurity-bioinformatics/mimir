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
    "translations",
    "hmmer_output",
    "hmm_max_evalue", 
    "hmm_min_score", 
    "hmm_max_hits",   
    "hmm_min_acc",   
    "hmm_max_gap"    
)
lapply(nf_vars, nf_var_check)

### process variables 
# read in fasta file as DNAbin
seqs <- ape::read.FASTA(fasta_file)

# read in translated sequences
translations <- ape::read.FASTA(translations, type = "AA")

# read in hmmer_output
domtblout <- readr::read_lines(hmmer_output, lazy=FALSE, progress=FALSE) 

# filtering parameters
evalue_threshold <- hmm_max_evalue %>% as.numeric()
score_threshold <- hmm_min_score %>% as.numeric()
domain_threshold <- hmm_max_hits %>% as.integer() # max allowed domains (1 is safe)
acc_threshold <- hmm_min_acc %>% as.numeric()
terminalgap_threshold <- hmm_max_gap %>% as.integer() # terminal gap is the distance between the beginning or end of the envelope/HMM and the min/max of the target or HMM
## NOTE: terminal gaps above the threshold on the same side (start or end) for both the envelope and HMM is expected in the case of frameshifts
## More extensive frameshift detection could be done with alignments to other sequences and removing those with gaps not in multiples of three

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
    # read in element as if it were a tsv file
    readr::read_tsv(
        col_names = names(col_types$cols),
        na = '-',
        col_types = col_types,
        lazy = FALSE, 
        progress = FALSE
    ) %>% 
    # replace '\t' in description column with spaces to regenerate sequence names
    dplyr::mutate(
        description = stringr::str_replace_all(description, "\\\t", " ")
    )

# convert to single sequence hits
hits <- 
    domtblout_parsed %>%
    # undo HMMER splitting sequence names into two parts at the first space character
    dplyr::mutate(
        target_name = dplyr::if_else(
            description != "-",
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
    # work out if any stop codons lie within the hit envelope
    tidyr::separate_longer_delim(
        cols = stop_codons,
        delim = "|"
    ) %>%
    dplyr::mutate(stop_codons = as.integer(stop_codons)) %>%
    dplyr::mutate(
        within_hit = dplyr::case_when(
            between(stop_codons, env_from, env_to) ~ TRUE,
            .default = FALSE
        )
    ) %>%
    dplyr::group_by(target_name, frame) %>%
    dplyr::mutate(
        stop_codon_hit = any(within_hit), # TRUE if any stop codons found in hit envelope
        stop_codons = paste(stop_codons, collapse = "|"), # "un-separate" the stop_codons column
        stop_codons = dplyr::na_if(stop_codons, "NA")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-within_hit) %>%
    dplyr::distinct() 

# check hits are all in sequence file
if(!all(hits$target_name %in% names(seqs))){
    stop("One or more sequence names in HMMER output table do not match names in input .fasta file")
}

# filter hits based on various thresholds
hits_retained <- 
    hits %>%
    # remove sequences that don't meet all criteria
    dplyr::filter(
        evalue_i <= evalue_threshold & 
        score >= score_threshold & 
        # domain_total <= domain_threshold & 
        acc >= acc_threshold &
        stop_codon_hit == FALSE &
        !( hmm_start_gap > terminalgap_threshold & env_start_gap > terminalgap_threshold ) & # both need to be fulfilled 
        !( hmm_end_gap > terminalgap_threshold & env_end_gap > terminalgap_threshold ) 
    )

if( nrow(hits_retained) > 0){

    # forward strand hits
    hits_fwd <- 
        hits_retained %>%
        dplyr::filter(frame > 0)

    # reverse strand hits
    hits_rev <- 
        hits_retained %>%
        dplyr::filter(frame < 0)

    # fwd seqs
    if (any(names(seqs) %in% hits_fwd$target_name)){
        seqs_fwd <- seqs[names(seqs) %in% hits_fwd$target_name]
    } else {
        seqs_fwd <- list() %>% as.DNAbin()
    }

    # rev seqs
    if (any(names(seqs) %in% hits_rev$target_name)){
        seqs_rev <- seqs[names(seqs) %in% hits_rev$target_name]
    } else {
        seqs_rev <- list() %>% as.DNAbin()
    }

    # revcomp the seqs with reverse strand hits
    seqs_revcomp <- 
        DNAbin2DNAstringset(seqs_rev) %>%
        Biostrings::reverseComplement(.) %>%
        ape::as.DNAbin(.)

    # combine seqs both in same orientation
    seqs_combined <- 
        concat_DNAbin(seqs_fwd, seqs_revcomp)

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

    ## subset translations based on location of envelope for each sequence (for further trimming with trimmed HMM)

    # remove translations not present in the subset sequences
    translations_retained <- translations[names(translations) %in% hit_locations$old_name]

    # get order of translations
    translations_retained_names <- names(translations_retained)

    # get padded envelope locations
    env_locations <- 
        hit_locations %>% 
        # pad envelope 1 residue either side
        dplyr::mutate(
            env_from_pad = base::pmax(1, env_from - 1),
            env_to_pad = base::pmin(target_len, env_to + 1)
        ) %>%
        # force order of tibble to be the same as the AAbin object
        dplyr::mutate( old_name = forcats::fct_relevel(old_name, translations_retained_names) ) %>% 
        dplyr::arrange(old_name)

    # subset translations to padded envelopes
    translations_subset <-
        translations_retained %>%
        # convert to AAstringset
        as.list() %>%
        as.character() %>%
        purrr::map(function(y){ paste0(y, collapse = "") }) %>%
        unlist %>%
        Biostrings::AAStringSet() %>%
        # subset
        XVector::subseq(., start = env_locations$env_from_pad, end = env_locations$env_to_pad) %>%
        ape::as.AAbin()

    ##
    # seqs without a significant hit (either no hit or an excluded hit)
    seqs_nohit <- seqs[!names(seqs) %in% c(hits_fwd$target_name, hits_rev$target_name)]

    # save data for hits that were removed after filtering
    seqs_removed_tibble <- 
        hits %>%
        # remove sequences that don't meet all criteria
        dplyr::filter(
            evalue_i > evalue_threshold |
            score < score_threshold | 
            # domain_total > domain_threshold |
            acc < acc_threshold |
            stop_codon_hit == TRUE |
            ( hmm_start_gap > terminalgap_threshold & env_start_gap > terminalgap_threshold ) | # both need to be fulfilled 
            ( hmm_end_gap > terminalgap_threshold & env_end_gap > terminalgap_threshold ) 
        ) %>%
        # which filters were failed
        dplyr::mutate(
            fail_max_evalue =       evalue_i > evalue_threshold,
            fail_min_score =        score < score_threshold,
            # fail_max_hits =         domain_total > domain_threshold ,
            fail_min_acc =          acc < acc_threshold ,
            fail_max_gap_start =    ( hmm_start_gap > terminalgap_threshold & env_start_gap > terminalgap_threshold ),
            fail_max_gap_end =      ( hmm_end_gap > terminalgap_threshold & env_end_gap > terminalgap_threshold ),
            fail_stop_codon =       stop_codon_hit == TRUE
        ) %>%
        dplyr::mutate(removal_type = "excluded_hit", .after = stop_codon_hit) %>%
        # add names of sequences without a hit
        tibble::add_row(
            target_name = names(seqs_nohit)[!names(seqs_nohit) %in% .$target_name], 
            removal_type = "no_hit"
        ) %>%
        dplyr::select(-old_name) # remove old_name as redundant with target_name

    readr::write_csv(seqs_removed_tibble, file = "removed_full.csv")

    ### outputs


    # check all sequences are either retained or removed
    if ( length(seqs_subset) + length(seqs_nohit) != length(seqs) ){
        stop("ERROR: Sum of retained and removed sequences does not equal the number of input sequences")
    }

    # write fasta of subsetted sequences
    if ( !is.null(seqs_subset) && length(seqs_subset) > 0 ){
        write_fasta(
            seqs_subset, 
            file = paste0("retained_full.fasta"), 
            compress = FALSE
            )
    } else {
        file.create(paste0("retained_full.fasta"))
    }

    # write fasta of excluded (filtered-out) sequences
    if ( !is.null(seqs_nohit) && length(seqs_nohit) > 0 ){
        write_fasta(
            seqs_nohit, 
            file = paste0("removed_full.fasta"), 
            compress = FALSE
            )
    } else {
        file.create(paste0("removed_full.fasta"))
    }

    # write trimmed translations
    if (!is.null(translations_subset) && length(translations_subset) > 0){
        ape::write.FASTA(
        translations_subset,
        file = paste0("translations_retained.fasta")
        )  
    } else {
        file.create(paste0("translations_retained.fasta"))
    }

} else {

    # seqs without a significant hit (either no hit or an excluded hit)
    seqs_nohit <- seqs

    # save data for hits that were removed after filtering
    seqs_removed_tibble <- 
        hits %>%
        # remove sequences that don't meet all criteria
        dplyr::filter(
            evalue_i > evalue_threshold |
            score < score_threshold | 
            # domain_total > domain_threshold |
            acc < acc_threshold |
            stop_codon_hit == TRUE |
            ( hmm_start_gap > terminalgap_threshold & env_start_gap > terminalgap_threshold ) | # both need to be fulfilled 
            ( hmm_end_gap > terminalgap_threshold & env_end_gap > terminalgap_threshold ) 
        ) %>%
        # which filters were failed
        dplyr::mutate(
            fail_max_evalue =       evalue_i > evalue_threshold,
            fail_min_score =        score < score_threshold,
            # fail_max_hits =         domain_total > domain_threshold ,
            fail_min_acc =          acc < acc_threshold ,
            fail_max_gap_start =    ( hmm_start_gap > terminalgap_threshold & env_start_gap > terminalgap_threshold ),
            fail_max_gap_end =      ( hmm_end_gap > terminalgap_threshold & env_end_gap > terminalgap_threshold ),
            fail_stop_codon =       stop_codon_hit == TRUE
        ) %>%
        dplyr::mutate(removal_type = "excluded_hit", .after = stop_codon_hit) %>%
        # add names of sequences without a hit
        tibble::add_row(
            target_name = names(seqs_nohit)[!names(seqs_nohit) %in% .$target_name], 
            removal_type = "no_hit"
        ) %>%
        dplyr::select(-old_name) # remove old_name as redundant with target_name

    readr::write_csv(seqs_removed_tibble, file = "removed_full.csv")

    write_fasta(
        seqs_nohit, 
        file = paste0("removed_full.fasta"), 
        compress = FALSE
    )

    file.create(paste0("retained_full.fasta"))
    file.create(paste0("translations_retained.fasta"))


}
