#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "Biostrings",
    "DECIPHER",
    "IRanges",
    "R.utils",
    "RCurl",
    "ape",
    # "aphid",
    "bold",
    "data.table",
    "data.tree",
    "dplyr",
    # "entropy",
    # "fs",
    # "furrr",
    # "future",
    "httr",
    "kmer",
    "magrittr",
    "methods",
    "openssl",
    "parallel",
    # "phytools",
    "purrr",
    "readr",
    # "rentrez",
    # "rvest",
    "stats",
    "stringr",
    # "taxize",
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
    "primer_fwd",
    "primer_rev",
    "remove_primers",
    "max_primer_mismatches",
    "min_length_trimmed"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# sequences aligned to primers using "mafft --add-fragments"
dap_alignment <- ape::read.FASTA(fasta_file) %>% DNAbin2DNAstringset(.)

# parameter 'remove_primers'
if (remove_primers == "true"){
    remove_primers <- TRUE
} else if (remove_primers == "false"){
    remove_primers <- FALSE
} else {
    stop(paste0("'remove_primers' value of '",remove_primers,"' should be 'true' or 'false'"))
}

# parameter 'max_primer_mismatches' -- maximum allowed mismatches to consensus of alignment for either of the best primers
# high values may give you incorrect results
max_primer_mismatches <- max_primer_mismatches %>% as.integer()

# minimum length of sequences to retain after trimming to primer region
min_length_trimmed <- min_length_trimmed %>% as.integer()

# minimum number of disambiguated primer sequences that need to support the location of a primer
location_threshold <- 0.8

# new ambiguity map that handles gaps 
NEW_AMB_MAP <- c(Biostrings::IUPAC_CODE_MAP, "-" = "-")

# get a logical vector of matches between the primer string and the consensus string, using bytes
# code from https://support.bioconductor.org/p/9159133/
# each element of vector represents a position on the string and is TRUE if IUPAC codes are congruent
BPComparison <- function(dna1, dna2){
    (as.raw(dna1) & as.raw(dna2))  > 0
}

BPComparisonNoGaps <- function(dna1, dna2){
    ((as.raw(dna1) & as.raw(dna2)) & as.raw(15))  > 0
}

### run code

# NOTE: 'primers_same' is fwd + rev_rc; 'primers_oppo' is rev + fwd_rc

# create DNAString objects for all four primer strings
fwd_ag.string <- primer_fwd
rev_ag.string <- primer_rev
fwd_ag.seq <- Biostrings::DNAString(fwd_ag.string)
rev_ag.seq <- Biostrings::DNAString(rev_ag.string)
fwd_rc.seq <- Biostrings::reverseComplement(fwd_ag.seq)
rev_rc.seq <- Biostrings::reverseComplement(rev_ag.seq)
fwd_rc.string <- fwd_rc.seq %>% as.character()
rev_rc.string <- rev_rc.seq %>% as.character()

strings_vec <- c(fwd_ag.string, rev_ag.string, fwd_rc.string, rev_rc.string)

# check primer sequences are all different
if (any(base::duplicated(strings_vec))){
    stop("Not all primer sequences are unique (comparing 'fwd', 'rev' and their reverse complements)")
}



# # nested list of primer sequences
# pairs_list <- list(
#     list(primers_same.seq, fwd_ag.string, rev_rc.string),
#     list(primers_oppo.seq, rev_ag.string, fwd_rc.string)
# )

# names(pairs_list) <- c("primers_same", "primers_oppo")

# # get a logical vector of matches between the primer string and the consensus string, using bytes
# # code from https://support.bioconductor.org/p/9159133/
# # each element of vector represents a position on the string and is TRUE if IUPAC codes are congruent
# BPComparison <- function(dna1, dna2){
#     (as.raw(dna1) & as.raw(dna2))  > 0
# }

# BPComparisonNoGaps <- function(dna1, dna2){
#     ((as.raw(dna1) & as.raw(dna2)) & as.raw(15))  > 0
# }

# # new ambiguity map that handles gaps 
# NEW_AMB_MAP <- c(IUPAC_CODE_MAP, "-" = "-")

# ## function to give consensus scores for each primer in a list
# primer_score <- function(x, idx, msa){
#     # 1-seq gapped alignment just of primer pair
#     pa <- msa[names(msa) == idx][[1]]
    
#     ## fwd primer
#     # pattern to match 1st primer with arbitrary gaps
#     p1_pattern <- 
#         x[[2]] %>% 
#         stringr::str_replace_all(., "(\\w)", "\\1(-)*") %>%
#         stringr::str_remove(., "\\(-\\)\\*$")
    
#     # location of 1st primer
#     p1_location <- 
#         pa %>%
#         as.character() %>%
#         stringr::str_locate(., p1_pattern) %>% .[1:2]
    
#     # subset database msa to primer positions and get consensus of non-primer sequences
#     p1.consensus <- 
#         msa %>% 
#         XVector::subseq(., start = p1_location[1], end = p1_location[2]) %>%
#         # remove primer sequences at the end
#         .[!names(msa) %in% c("primers_same", "primers_oppo")] %>%
#         # get consensus of region
#         Biostrings::consensusMatrix(., as.prob = FALSE) %>%
#         Biostrings::consensusString(ambiguityMap = NEW_AMB_MAP, threshold = 0.1) %>%
#         Biostrings::DNAString()
    
#     # get sum of matches between the primer string and the consensus string, using bytes
#     # gaps are treated as non-matches
#     p1.bpcomp <- 
#         BPComparisonNoGaps(Biostrings::DNAString(x[[2]]), p1.consensus) %>%
#         # sum up the number of correct matches
#         sum() 
    
#     # no. of mismatches to primer string
#     p1.bpmm <- width(x[[2]]) - p1.bpcomp
    
#     ## 2nd primer
#     # pattern to match 2nd primer with arbitrary gaps
#     p2_pattern <- 
#         x[[3]] %>% 
#         stringr::str_replace_all(., "(\\w)", "\\1(-)*") %>%
#         stringr::str_remove(., "\\(-\\)\\*$")
    
#     # location of 2nd primer
#     p2_location <- 
#         pa %>%
#         as.character() %>%
#         stringr::str_locate(., p2_pattern) %>% .[1:2]
    
#     # subset database msa to primer positions and get consensus of non-primer sequences
#     p2.consensus <- 
#         msa %>% 
#         XVector::subseq(., start = p2_location[1], end = p2_location[2]) %>%
#         # remove primer sequences at the end
#         .[!names(msa) %in% c("primers_same", "primers_oppo")] %>%
#         # get consensus of region
#         Biostrings::consensusMatrix(., as.prob = FALSE) %>%
#         Biostrings::consensusString(ambiguityMap = NEW_AMB_MAP, threshold = 0.1) %>%
#         Biostrings::DNAString()
    
#     # get sum of matches between the primer string and the consensus string, using bytes
#     # gaps are treated as non-matches
#     p2.bpcomp <- 
#         BPComparisonNoGaps(Biostrings::DNAString(x[[3]]), p2.consensus) %>%
#         # sum up the number of correct matches
#         sum() 
    
#     # no. of mismatches to primer string
#     p2.bpmm <- width(x[[3]]) - p2.bpcomp
    
#     return(c(p1.bpmm, p2.bpmm))
# }

# # get primer mismatch scores
# pms <- 
#     purrr::imap(
#         pairs_list,
#         primer_score,
#         msa = seqs_dss
#     ) 

# # calculate combined score for pairs
# scores_same <- pms$primers_same
# scores_oppo <- pms$primers_oppo

# print(paste0("Same orientation scores: ", stringr::str_c(scores_same, collapse = " and ")))
# print(paste0("Opposite orientation scores: ", stringr::str_c(scores_oppo, collapse = " and ")))

# # work out best scoring pair
# if ( sum(scores_same) < sum(scores_oppo) ){
#     bsp <- "same"
#     bsp_scores <- scores_same
# } else if ( sum(scores_same) > sum(scores_oppo) ){
#     bsp <- "oppo"
#     bsp_scores <- scores_oppo
# } else {
#     stop("Combined mismatch scores for both primer orientations are the same; cannot choose best primer orientation!")
# }

# # check neither best scoring primer has more mismatches to consensus than allowed
# if ( any(bsp_scores) > max_primer_mismatches ){
#     print(bsp_scores)
#     stop(paste0("At least one mismatch score for the best scoring primer orientation is greater than '--max_primer_mismatch'"))
# }

# # remove losing primer pair, rename winning pair, change orientation if needed, and remove columns that only contain gaps
# # if winning primer pair is "oppo", reverse complement the aligned database so the fwd primer is upstream/left
# if ( bsp == "same"){
#     seqs_dss_chosen <- 
#         seqs_dss %>%
#         .[!names(.) == "primers_oppo"] %>%
#         DECIPHER::RemoveGaps(., removeGaps = "common")
#     # rename primers
#     names(seqs_dss_chosen) <- names(seqs_dss_chosen) %>% replace(., . == "primers_same", "primers")
# } else if ( bsp == "oppo" ){
#     seqs_dss_chosen <- 
#         Biostrings::reverseComplement(seqs_dss) %>%
#         .[!names(.) == "primers_same"] %>%
#         DECIPHER::RemoveGaps(., removeGaps = "common")
#     # rename primers
#     names(seqs_dss_chosen) <- names(seqs_dss_chosen) %>% replace(., . == "primers_oppo", "primers")
# } else {
#     stop(paste0("'bsp' value of '",bsp,"' is unexpected"))
# }
# rm(seqs_dss)
# gc()

# ### trim database sequences to primer positions
# # because if primers_oppo won the database was revcomped, primer sequences are fwd_ag and rev_rc no matter what

# # alignment of primers
# primer.alignment <- seqs_dss_chosen[names(seqs_dss_chosen) == "primers"][[1]]

# ## position of fwd primer
# # pattern to match fwd primer with arbitrary gaps
# fwd_pattern <- 
#     fwd_ag.string %>% 
#     stringr::str_replace_all(., "(\\w)", "\\1(-)*") %>%
#     stringr::str_remove(., "\\(-\\)\\*$")

# # get location of fwd
# fwd_location <- 
#     primer.alignment %>%
#     as.character() %>%
#     stringr::str_locate(., fwd_pattern) %>% .[1:2]

# ## position of rev primer
# # pattern to match rev primer with arbitrary gaps
# rev_pattern <- 
#     rev_rc.string %>% 
#     stringr::str_replace_all(., "(\\w)", "\\1(-)*") %>%
#     stringr::str_remove(., "\\(-\\)\\*$")

# # get location of fwd
# rev_location <- 
#     primer.alignment %>%
#     as.character() %>%
#     stringr::str_locate(., rev_pattern) %>% .[1:2]

# if (remove_primers) {
#     # included region does not include primers
#     trim_from <- fwd_location[2] + 1
#     trim_to <- rev_location[1] - 1
# } else {
#     # included region does include primers
#     trim_from <- fwd_location[1]
#     trim_to <- rev_location[2]
# }

# # subset alignment by positions
# seqs_trimmed <- 
#     seqs_dss_chosen %>%
#     XVector::subseq(., start = trim_from, end = trim_to) %>%
#     # remove primers from alignment
#     .[!names(.) == "primers"]

# # filter out sequences shorter than '--min_length_trimmed'
# # convert to ungapped, and get names of sequences shorter than threshold
# below_min_length <- 
#     seqs_trimmed %>%
#     DECIPHER::RemoveGaps(., removeGaps = "all") %>%
#     width(.) < min_length_trimmed

# # remove short sequences
# seqs_pass_length <- 
#     seqs_trimmed[!below_min_length]

# # write sequences out
# ape::write.FASTA(seqs_pass_length %>% ape::as.DNAbin(), file = "trimmed.fasta")

# ## write out sequences that fail length after trimming
# ape::write.FASTA(seqs_trimmed[below_min_length] %>% ape::as.DNAbin(), file = "trimmed.removed.fasta")








#### new consensus implementation


## getting location of degenerate primer, then # mismatches to alignment consensus
# nested list of primer name and sequence string
primers_list <- 
    list(
        list("fwd_ag", fwd_ag.string),
        list("rev_ag", rev_ag.string),
        list("fwd_rc", fwd_rc.string),
        list("rev_rc", rev_rc.string)
    )

# function returns list of:
# [[1]] primer name
# [[2]] proportion of DA sequences at most common position (double)
# [[3]] most common position (vector of start and end)
# [[4]] number of mismatches to the alignment consensus at that position (integer)
# primer = element of primer_list
# msa = complete alignment including primers 
primer_check <- function(primer, msa, location_threshold = location_threshold){
  
    # DA sequence alignments for given primer
    p.dap <- msa[names(msa) %>% stringr::str_starts(primer[[1]])]
    
    # location of each DA sequence
    p.locs <- 
    lapply(
        p.dap, # apply to each sequence of DSS object
        function(x){
            # find gaps
            da.gap_pos <- 
                x %>%
                as.character() %>%
                stringr::str_locate_all(., "-+") %>% .[[1]]
            # find start and end of nucleotides
            da.start <- da.gap_pos[1,2] + 1
            da.start <- unname(da.start)
            da.end <- da.gap_pos[nrow(da.gap_pos),1] - 1
            da.end <- unname(da.end)
            # return positions separated by /
            return(paste0(da.start, "/", da.end))
        }
    )
    
    # determine which position wins
    p.position <-
        p.locs %>% 
        unlist() %>% 
        table() %>% # get counts of each position
        tibble::as_tibble() %>%
        dplyr::rename(., location = `.`) %>%
        dplyr::mutate(
            prop = n / sum(n) # proportion of DA sequences that support this location
        ) %>%
        dplyr::arrange(desc(prop)) %>% # arrange by highest prop
        dplyr::slice(1) # get highest rating position

    # print detail for log
    print(
        paste0(
            "Primer '",
            primer[[1]],
            "': ", 
            round(p.position$prop * 100, 1), 
            "% (", 
            p.position$n, 
            ") sequences agree on position '", 
            p.position$location,
            "'"
        )
    )
    
    # get start and end positions as integer vector
    p.posvec <-
        p.position %>% 
        dplyr::pull(location) %>%
        stringr::str_split(., "/", n = 2) %>%
        unlist() %>% 
        as.integer()
    
    # check winning position is a vector of two values (start and end)
    if (p.posvec %>% length() != 2){
        stop(paste0("Could not get consensus position for '",primer[[1]],"' primer"))
    }
    
    # subset database to positions
    p.cons <- 
        msa %>%
        XVector::subseq(., start = p.posvec[1], end = p.posvec[2]) %>%
        # remove primer sequences at the end
        .[!names(msa) %>% stringr::str_starts(., "fwd_ag|rev_ag|fwd_rc|rev_rc")] %>%
        # get consensus of region
        Biostrings::consensusMatrix(., as.prob = FALSE) %>%
        Biostrings::consensusString(ambiguityMap = NEW_AMB_MAP, threshold = 0.1) %>%
        Biostrings::DNAString()
    
    # get sum of matches between the primer string and the consensus string, using bytes
    # gaps are treated as non-matches
    bpcomp <- 
        BPComparisonNoGaps(Biostrings::DNAString(primer[[2]]), p.cons) %>%
        # sum up the number of correct matches
        sum() 
    
    # no. of mismatches to primer string
    bpmm <- width(primer[[2]]) - bpcomp
    
    # return vector
    return(
        c(
            primer[[1]], # primer name
            p.position$prop, # prop sequences supporting position
            p.posvec, # vector of start and end positions
            bpmm # number of mismatches to alignment consensus
        )
    )
}

# apply function and format output into tibble
primer_check_output <-
    lapply(
        primers_list, 
        primer_check, 
        msa = dap_alignment
    ) %>% 
    as_tibble(.name_repair = "unique") %>%
    dplyr::mutate(
        type = c("name","loc_agree", "start", "end", "mismatches")
    ) %>%
    tidyr::pivot_longer(...1:...4, names_to = "old") %>%
    tidyr::pivot_wider(names_from = type, values_from = value) %>%
    dplyr::select(-old) %>%
    dplyr::mutate(
        loc_agree = as.double(loc_agree),
        start = as.integer(start),
        end = as.integer(end),
        mismatches = as.integer(mismatches)
    )

# get mismatch scores for each primer pair
scores_same <- c(
    primer_check_output %>% dplyr::filter(name == "fwd_ag") %>% dplyr::pull(mismatches),
    primer_check_output %>% dplyr::filter(name == "rev_rc") %>% dplyr::pull(mismatches)
)
scores_oppo <- c(
    primer_check_output %>% dplyr::filter(name == "fwd_rc") %>% dplyr::pull(mismatches),
    primer_check_output %>% dplyr::filter(name == "rev_ag") %>% dplyr::pull(mismatches)
)

print(paste0("Same orientation scores: ", stringr::str_c(scores_same, collapse = " and ")))
print(paste0("Opposite orientation scores: ", stringr::str_c(scores_oppo, collapse = " and ")))

# work out best scoring pair
if ( sum(scores_same) < sum(scores_oppo) ){
    bsp <- "same"
    bsp_scores <- scores_same
    pco_best <- 
        primer_check_output %>%
        dplyr::filter(name %in% c("fwd_ag","rev_rc")) %>%
        dplyr::mutate(
            name = dplyr::case_when(
                name == "fwd_ag" ~ "fwd",
                name == "rev_rc" ~ "rev"
            )
        )
} else if ( sum(scores_same) > sum(scores_oppo) ){
    bsp <- "oppo"
    bsp_scores <- scores_oppo
    pco_best <- 
        primer_check_output %>%
        dplyr::filter(name %in% c("fwd_rc","rev_ag")) %>%
        dplyr::mutate(
            name = dplyr::case_when(
                name == "fwd_rc" ~ "fwd",
                name == "rev_ag" ~ "rev"
            )
        )
} else {
    stop("Combined mismatch scores for both primer orientations are the same; cannot choose best primer orientation!")
}

if (!bsp %in% c("same","oppo")){
    stop("bsp not valid")
}

# check best scoring pair has confident locations
if ( any(pco_best$loc_agree < location_threshold) ){
    stop(
        paste0(
            "Best scoring primer pair has consensus locations of ", 
            stringr::str_flatten(round(pco_best$loc_agree, 2), collapse = " and "),
            ", when minimum threshold is ",
            location_threshold,
            ". \nPrimer sequences may be invalid or incorrect"
        )
    )
}

# check best scoring pair has allowed number of mismatches to consensus
if ( any(bsp_scores) > max_primer_mismatches ){
    print(paste0("Scores: ",stringr::str_flatten(bsp_scores, collapse = " and ")))
    stop(paste0("At least one mismatch score for the best scoring primer orientation is greater than '--max_primer_mismatch'"))
}

# for alignment, remove primers and remove columns that only contain gaps (ie. were induced by primer alignment)
np_alignment <- 
    dap_alignment %>%
    .[!names(.) %>% stringr::str_starts(., "fwd_ag|rev_ag|fwd_rc|rev_rc")] %>%
    DECIPHER::RemoveGaps(., removeGaps = "common")

### trim database sequences to primer positions
# use existing location info 
if (remove_primers) {
    # included region does not include primers
    trim_from <- pco_best[pco_best$name == "fwd",]$end + 1
    trim_to <- pco_best[pco_best$name == "rev",]$start - 1
} else {
    # included region does include primers
    trim_from <- pco_best[pco_best$name == "fwd",]$start
    trim_to <- pco_best[pco_best$name == "rev",]$end
}

# subset alignment by positions
alignment_trimmed <- 
  np_alignment %>%
  XVector::subseq(., start = trim_from, end = trim_to)

# reverse orientation of database if "oppo" was best orientation of primers
if (bsp == "oppo"){
    alignment_trimmed <- alignment_trimmed %>% Biostrings::reverseComplement()
} 

# filter out sequences shorter than '--min_length_trimmed'
# convert to ungapped, and get names of sequences shorter than threshold
below_min_length <- 
    alignment_trimmed %>%
    DECIPHER::RemoveGaps(., removeGaps = "all") %>%
    width(.) < min_length_trimmed

# remove short sequences
alignment_final <- alignment_trimmed[!below_min_length]

# write trimmed sequences out
ape::write.FASTA(alignment_final %>% ape::as.DNAbin(), file = "trimmed.fasta")

## write out sequences that fail length after trimming
ape::write.FASTA(alignment_trimmed[below_min_length] %>% ape::as.DNAbin(), file = "trimmed.removed.fasta")

