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
    "fasta_file"
)
lapply(nf_vars, nf_var_check)

### process variables 

# read in fasta file as DNAbin
seed_primers <- ape::read.FASTA(fasta_file, type = "AA")

# convert to AAStringSet
seed_primers.aass <- 
    seed_primers %>%
    as.list() %>%
    as.character() %>%
    purrr::map(function(y){ paste0(y, collapse = "") }) %>%
    unlist %>%
    Biostrings::AAStringSet()

# read in primers file as DNAbin
primers_seq <- ape::read.FASTA(primers_file, type = "DNA")

# convert to DNAStringSet
primers_seq.dss <- 
    primers_seq %>%
    as.list() %>%
    as.character() %>%
    purrr::map(function(y){ paste0(y, collapse = "") }) %>%
    unlist %>%
    Biostrings::DNAStringSet()

### run code

# four primer orientations, fwd_ag, rev_ag, fwd_rc, rev_rc
# three frames per primer: 1, 2, 3
# translated sequences have been aligned without adding additional gaps, so may be missing some sequence

primers <- c("fwd_ag","rev_ag","fwd_rc","rev_rc")

# just seeds without primers
seed.aass <- seed_primers.aass[!names(seed_primers.aass) %>% stringr::str_detect(., "_frame=\\d$")]

# just primers without seeds
primers.aass <- seed_primers.aass[names(seed_primers.aass) %>% stringr::str_detect(., "_frame=\\d$")]

# consensus matrix of seeds
seed.cm <- Biostrings::consensusMatrix(seed.aass, as.prob = T) 

# per frame across all primers
frame_agreement <- 
    lapply(
        primers,
        function(x){
            lapply(
                c("1","2","3"), 
                function(frame, primer = x, msa = primers.aass){
                    # just alignments of primers for primer x frame combo
                    pf_align <- msa[names(msa) %>% stringr::str_detect(., paste0("^",primer,"\\d+_frame=",frame,"$"))]
                    # per alignment in primer x frame combo
                    lapply(
                        1:length(pf_align),
                        function(i){
                            palign <- pf_align[[i]]
                            pname <- names(pf_align[i])
                            prep <- stringr::str_extract(pname, paste0("(?<=", primer, ")\\d+(?=_frame\\=",frame,")"))
                            # positions of alignment characters
                            ppos <- 
                                palign %>%
                                as.character() %>%
                                stringr::str_locate_all(., "[^-]") %>% .[[1]] %>% .[,1]
                            pstart <- min(ppos)
                            pend <- max(ppos)
                            # get mean probability of primer alignment
                            pmean <- 
                                # get probability for each position
                                lapply(
                                    ppos, 
                                    function(x, cm = seed.cm, pa = palign){
                                        # primer character at position
                                        pchar <- as.character(pa[x])
                                        # probability of primer character at position in seed alignment
                                        aprob <- cm[rownames(cm) == pchar,x] %>% unname()
                                        return(aprob)
                                    }
                                ) %>% 
                                unlist() %>% 
                                # get mean probability across all positions
                                mean()
                            # output tibble of info per non-degenerate sequence
                            out <- tibble::tibble(
                                primer = primer, 
                                rep = prep, 
                                frame = frame,
                                start = pstart,
                                end = pend,
                                location = paste0(start,"-",end),
                                agreement = pmean
                            )
                            return(out)
                        }
                    ) %>%
                        dplyr::bind_rows() %>%
                        dplyr::summarise(
                            .by = c(primer, frame, location, start, end),
                            agreement = mean(agreement),
                            n = n()
                        ) %>%
                        dplyr::mutate(prop = n/sum(n)) %>%
                        # sort by plurality location, break ties with agreement
                        dplyr::arrange(desc(n), desc(agreement)) %>%
                        dplyr::slice(1) 
                }
            ) %>% 
                dplyr::bind_rows()
        }
    ) %>% 
    dplyr::bind_rows()

# best frame per primer (based on agreement; prop as tiebreaker)
best_frames <- 
    frame_agreement %>%
    dplyr::group_by(primer) %>%
    dplyr::arrange(desc(agreement), desc(prop), .by_group = T) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

# get mean agreement for each primer orientation
agree_same <- best_frames %>% dplyr::filter(primer %in% c("fwd_ag", "rev_rc")) %>% dplyr::pull(agreement) %>% mean
agree_oppo <- best_frames %>% dplyr::filter(primer %in% c("fwd_rc", "rev_ag")) %>% dplyr::pull(agreement) %>% mean

# choose best primer orientation 
if ( agree_same > agree_oppo ){
    n_start <- "fwd_ag"
    n_end <- "rev_rc"
} else if ( agree_same < agree_oppo ){
    n_start <- "rev_ag"
    n_end <- "fwd_rc"
} else {
    stop("Scores for each primer orientation are the same, unable to decide")
}

# and determine start and end positions in alignment, plus primer offsets and primer lengths
pos_start <- best_frames %>% dplyr::filter(primer == n_start) %>% dplyr::pull(start)
pos_end <- best_frames %>% dplyr::filter(primer == n_end) %>% dplyr::pull(end)

frame_start <- best_frames %>% dplyr::filter(primer == n_start) %>% dplyr::pull(frame) %>% as.numeric() 
plen_start <- primers_seq.dss[names(primers_seq.dss) == n_start][[1]] %>% length()
# start offset is just the frame - 1
offset_start <- -(frame_start - 1)

frame_end <- best_frames %>% dplyr::filter(primer == n_end) %>% dplyr::pull(frame) %>% as.numeric() 
plen_end <- primers_seq.dss[names(primers_seq.dss) == n_end][[1]] %>% length()
# end offset is the remainder of the primer after the frame offset and the translated portion of the primer
offset_end <- (plen_end - (frame_end - 1)) %% 3

# check start position is smaller than end position
if ( pos_start > pos_end ){
    stop("Primer positions not possible, check code")
}

# output csv with frame offset and primer lengths
tibble::tibble(
    offset_start = offset_start,
    offset_end = offset_end,
    plen_start = plen_start,
    plen_end = plen_end
) %>%
    readr::write_csv(., "primer_info.csv")

# trim seed alignment to positions
seed_primers.aass %>%
    XVector::subseq(., start = pos_start, end = pos_end) %>%
    # remove primer sequences
    .[!names(.) %>% stringr::str_starts(., "fwd_ag|rev_ag|fwd_rc|rev_rc")] %>%
    ape::as.AAbin() %>%
    # write out trimmed alignment
    ape::write.FASTA(., "seed_trimmed.fasta")