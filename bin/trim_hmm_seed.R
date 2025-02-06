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

### run code

# four primer orientations, fwd_ag, rev_ag, fwd_rc, rev_rc
# three frames per primer: 1, 2, 3
# translated sequences have been aligned without adding additional gaps, so may be missing some sequence

primer_trans_check <- function(primer, msa, location_threshold = 0.5){
    frame <- c("1","2","3")
    # apply rest of function to each frame
    each_frame <- lapply(
            frame,
            function(frame, primer, msa){
                # sequence alignments for given primer and frame combo
                p.dap <- msa[names(msa) %>% stringr::str_detect(., paste0("^",primer[[1]],"\\d+_frame=",frame[[1]],"$"))]
                # location of each translated primer sequence
                p.locs <- 
                    lapply(
                        p.dap, # apply to each sequence of ASS object
                        function(x){
                            # find gaps
                            da.gap_pos <- 
                                x %>%
                                as.character() %>%
                                stringr::str_locate_all(., "-+") %>% .[[1]]
                            # find start and end of residues
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
                        prop = n / sum(n),
                        pass_threshold = prop >= {{location_threshold}}
                    ) %>%
                    dplyr::arrange(desc(prop)) %>%
                    dplyr::slice(1) # get highest rating position
                
                # print detail for log
                print(paste0(
                    "Primer '",
                    primer[[1]],
                    "', frame '",
                    frame[[1]],
                    ": ", 
                    round(p.position$prop * 100, 1), 
                    "% (", 
                    p.position$n, 
                    ") sequences agree on position '", 
                    p.position$location,
                    "'"
                ))
                
                # get start and end positions as integer vector
                p.posvec <-
                    p.position %>% 
                    dplyr::pull(location) %>%
                    stringr::str_split(., "/", n = 2) %>%
                    unlist() %>% 
                    as.integer()
                
                # check winning position is a vector of two values (start and end)
                if (p.posvec %>% length() != 2){
                    stop(paste0("Could not get consensus position for '",primer,"' primer"))
                }
                
                # return list
                return(
                    c(
                        primer,
                        frame,
                        p.position$prop,
                        p.posvec
                    )
                )
            }, 
            primer = primer, 
            msa = msa
    )
    frame_info <- 
        each_frame %>%
        tibble::as_tibble(.name_repair = "unique") %>%
        dplyr::mutate(
            type = c("name","frame", "prop", "start", "end")
        ) %>%
        tidyr::pivot_longer(...1:...3, names_to = "old") %>%
        tidyr::pivot_wider(names_from = type, values_from = value) %>%
        dplyr::select(-old) %>%
        dplyr::mutate(
            prop = as.double(prop),
            start = as.integer(start),
            end = as.integer(end),
            frame = as.integer(frame)
        ) 
    return(frame_info)
}

ptc_out <- 
    lapply(
        c("fwd_ag","rev_ag","fwd_rc","rev_rc"),
        primer_trans_check,
        msa = seed_primers.aass
    )

best_frames <-   
    ptc_out %>% 
    dplyr::bind_rows() %>%
    dplyr::group_by(name) %>%
    dplyr::arrange(desc(prop), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

prop_same <- best_frames %>% dplyr::filter(name %in% c("fwd_ag","rev_rc")) %>% dplyr::pull(prop) %>% mean(.)
prop_oppo <- best_frames %>% dplyr::filter(name %in% c("fwd_rc","rev_ag")) %>% dplyr::pull(prop) %>% mean(.)

if ( prop_same > prop_oppo ){
    best_orient <- best_frames %>% dplyr::filter(name %in% c("fwd_ag","rev_rc"))
    pos_start <- best_orient %>% dplyr::filter(name == "fwd_ag") %>% dplyr::pull(start)
    pos_end <- best_orient %>% dplyr::filter(name == "rev_rc") %>% dplyr::pull(end)
} else if ( prop_same < prop_oppo ){
    best_orient <- best_frames %>% dplyr::filter(name %in% c("fwd_rc","rev_ag"))
    pos_start <- best_orient %>% dplyr::filter(name == "rev_ag") %>% dplyr::pull(start)
    pos_end <- best_orient %>% dplyr::filter(name == "fwd_rc") %>% dplyr::pull(end)
} else {
    stop("Scores for each primer orientation are the same, unable to decide")
}

# check start position is smaller than end position
if ( pos_start > pos_end ){
    stop("Primer positions not possible, check code")
}

# trim seed alignment to positions
seed_primers.aass %>%
    XVector::subseq(., start = pos_start - 1, end = pos_end + 1) %>%
    # remove primer sequences
    .[!names(.) %>% stringr::str_starts(., "fwd_ag|rev_ag|fwd_rc|rev_rc")] %>%
    ape::as.AAbin() %>%
    # write out trimmed alignment
    ape::write.FASTA(., "seed_trimmed.fasta")