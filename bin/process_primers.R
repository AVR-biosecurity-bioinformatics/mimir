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
    "primer_fwd",
    "primer_rev",
    "gencodes_file",
    "marker_type"
    )
lapply(nf_vars, nf_var_check)

### process variables 

# import gencodes tibble
gencodes <- readr::read_delim(gencodes_file, delim = ",")

# code from modelr::typical 
find_typical <- function(x, ...) {
    counts <- table(x, exclude = NULL)
    names(counts)[max(counts) == counts]
}

# select genetic code type based on marker type and then the most frequent code 
genetic_code <- 
    gencodes %>%
    dplyr::mutate(
        code = dplyr::case_when(
            {{marker_type}} == "nuclear" ~ gencode,
            {{marker_type}} == "mitochondrial" ~ mtgencode,
            {{marker_type}} == "plastid" ~ plgencode,
            {{marker_type}} == "hydrogen" ~ hdgencode
        )
    ) %>%
    dplyr::pull(code) %>%
    as.character() %>%
    find_typical()

# check only one genetic code
if ( length(genetic_code) != 1 ){
    stop("Exactly one genetic code must be 'genetic_code' object")
}

### run code

# get sequences of all four primers
fwd_ag.string <- primer_fwd
rev_ag.string <- primer_rev
fwd_ag.seq <- Biostrings::DNAString(fwd_ag.string)
rev_ag.seq <- Biostrings::DNAString(rev_ag.string)
fwd_rc.seq <- Biostrings::reverseComplement(fwd_ag.seq)
rev_rc.seq <- Biostrings::reverseComplement(rev_ag.seq)

# save original primer sequences
primers.dss <- Biostrings::DNAStringSet(list(fwd_ag.seq, rev_ag.seq, fwd_rc.seq, rev_rc.seq))
names(primers.dss) <- c("fwd_ag","rev_ag", "fwd_rc", "rev_rc")
ape::as.DNAbin(primers.dss) %>% 
    ape::write.FASTA(., "primers_original.fasta")

## disambiguated sequences, named with degenerate primer name and a number from 1 to length()
# fwd_ag
fwd_ag.da <- fwd_ag.seq %>% DNAStringSet(.) %>% DECIPHER::Disambiguate(.) %>% .[[1]]
names(fwd_ag.da) <- seq(1,length(fwd_ag.da)) %>% stringr::str_replace(., "^", "fwd_ag")

# rev_ag
rev_ag.da <- rev_ag.seq %>% DNAStringSet(.) %>% DECIPHER::Disambiguate(.) %>% .[[1]]
names(rev_ag.da) <- seq(1,length(rev_ag.da)) %>% stringr::str_replace(., "^", "rev_ag")

# fwd_ag
fwd_rc.da <- fwd_rc.seq %>% DNAStringSet(.) %>% DECIPHER::Disambiguate(.) %>% .[[1]]
names(fwd_rc.da) <- seq(1,length(fwd_rc.da)) %>% stringr::str_replace(., "^", "fwd_rc")

# rev_rc
rev_rc.da <- rev_rc.seq %>% DNAStringSet(.) %>% DECIPHER::Disambiguate(.) %>% .[[1]]
names(rev_rc.da) <- seq(1,length(rev_rc.da)) %>% stringr::str_replace(., "^", "rev_rc")

# save primers as DSS
primers_disamb <- c(fwd_ag.da, rev_ag.da, fwd_rc.da, rev_rc.da)

# write to .fasta
primers_disamb %>% 
    ape::as.DNAbin() %>% 
    ape::write.FASTA(., "primers_disambiguated.fasta")

## translate primers (forward frames only)
frames <- c(
    lapply(1:3, function(pos) XVector::subseq(primers_disamb, start=pos))    )

# restructure list so each element is a DNAStringSet object of each sequence's 3 frames
frames_list <- 
    frames %>% 
    as.data.frame(col.names = c("f1","f2","f3")) %>% 
    tibble::as_tibble(rownames = "name") %>%
    tidyr::pivot_longer(
        cols = f1:f3, 
        names_to = "frame", 
        names_prefix = "f",
        values_to = "sequence"
    ) %>%
    dplyr::mutate(
        frame = stringr::str_replace(frame, "\\.","-") 
    ) %>%
    split(., .$name) %>% # split into one df per input sequence
    lapply(
        ., 
        function(x) {
            x_df <- as.data.frame(x)
            seqs_DSS <- DNAStringSet(x_df$sequence)
            names(seqs_DSS) <- x$name %>% stringr::str_replace(., "$", paste0("_frame=",x$frame)) # add frame to name
            return(seqs_DSS)
        } 
    )

# translate frames using genetic codes vector
translated <- 
    lapply( # applies function using combination of two vectors 
      frames_list,  
      function(y, z) {
            suppressWarnings(
                Biostrings::translate(
                    x = y, 
                    genetic.code = Biostrings::getGeneticCode(z),
                    if.fuzzy.codon=c("solve", "X"),
                    no.init.codon = TRUE
                )
            )
        }, 
        z = genetic_code
    )

# create a single object from the list of translated sequences
translations_combined <- AAStringSet()
for (i in 1:length(translated)){
    translations_combined <- S4Vectors::append(translations_combined, translated[[i]], after = length(translations_combined))
}

# convert stop codons to X character
translations_combined <- Biostrings::chartr("*", "X", translations_combined)

# save as .fasta
translations_combined %>% 
    ape::as.AAbin() %>%
    ape::write.FASTA(., file = "translated.fasta")
