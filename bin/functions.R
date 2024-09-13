#### custom R functions to be used in the pipeline

## checks existence of a variable and print error message if not defined
nf_var_check <- function(x) {
  if (!exists(x)) {
    stop(paste0("The variable'",x,"' is not defined! Make sure to check the Nextflow process inputs."))
  } else {
    print(paste0("Input variable '",x,"' = ",eval(parse(text = x))))
  }
}

## collapses repetive Groovy list variables down to a single variable
parse_nf_var_repeat <- function(x) {
  variable <- 
    stringr::str_extract_all(
      x, 
      pattern = "[^\\s,\\[\\]]+" # extract all runs of characters that aren't ' ' ',' '[' or ']' 
      ) %>% 
    unlist() %>%
    tibble::as_tibble_col(column_name = "col") %>% 
    unique() %>%
    dplyr::pull(col)
  
  if (length(variable) == 1) {
    out <- variable
  } else {
    out <- stop("*** nf variable contains multiple unique values! ***")
  }
  return(out)
}


#### taxreturn functions, ported over to allow modification
### 'depreciated' and 'data' functions not imported

# DNAbin Utilities -----------------------------------------------------

#' Convert DNABin to DNAStringSet
#'
#' @param x a DNABin object
#' @param remove_gaps Whether gaps should be removed
#'
#' @return
#' @export
#' @import purrr
#' @importFrom Biostrings DNA_ALPHABET
#' @importFrom Biostrings DNAStringSet
#' @importFrom methods is
#'
DNAbin2DNAstringset <- function (x, remove_gaps = FALSE) {
  if(!methods::is(x, "DNAbin")){
    stop("Input must be a DNAbin")
  }
  x %>%
    as.list() %>%
    as.character %>%
    purrr::map(function(y){
      y[!y %in% tolower(Biostrings::DNA_ALPHABET)] <- "N"
      if(isTRUE(remove_gaps)){
        y[y=="-"] <- ""
      }
      paste0(y, collapse="")
    })%>%
    unlist %>%
    Biostrings::DNAStringSet()
}

#' char2DNAbin
#'
#' @param x a character vector
#'
#' @return
#' @export
#'
#' @examples
char2DNAbin <- function (x) {
  dbytes <- as.raw(c(136, 24, 72, 40, 96, 144, 192, 48, 80,
                     160, 112, 224, 176, 208, 240, 240, 4, 2))
  indices <- c(65, 84, 71, 67, 83, 87, 82, 89, 75, 77, 66,
               86, 72, 68, 78, 73, 45, 63)
  vec <- raw(89)
  vec[indices] <- dbytes
  s2d1 <- function(s) vec[as.integer(charToRaw(s))]
  out <- lapply(x, s2d1)
  attr(out, "rerep.names") <- attr(x, "rerep.names")
  attr(out, "rerep.pointers") <- attr(x, "rerep.pointers")
  class(out) <- "DNAbin"
  return(out)
}

#' DNAbin to char
#'
#' @param x A DNAbin object
#'
#' @return
#' @export
#'
#' @examples
DNAbin2char <- function(x){
  cbytes <- as.raw(c(65, 84, 71, 67, 83, 87, 82, 89, 75, 77,
                     66, 86, 72, 68, 78, 45, 63))
  indices <- c(136, 24, 72, 40, 96, 144, 192, 48, 80, 160,
               112, 224, 176, 208, 240, 4, 2)
  vec <- raw(240)
  vec[indices] <- cbytes
  if (is.list(x)) {
    out <- sapply(x, function(s) rawToChar(vec[as.integer(s)]))
  } else {
    out <- rawToChar(vec[as.integer(x)])
  }
  attr(out, "rerep.names") <- attr(x, "rerep.names")
  attr(out, "rerep.pointers") <- attr(x, "rerep.pointers")
  return(out)
}

# add a queit arg
# add a detect "gz" set compress to true

#' write fasta
#'
#' @param x a list of sequences in DNAbin or AAbin format, or a vector of sequences as concatenated upper-case character strings.
#' @param file character string giving a valid file path to output the text to. If file = "" (default setting) the text file is written to the console.
#' @param compress logical indicating whether the output file should be gzipped.
#' @param quiet Whether progress should be printed to consoe
#'
#' @return
#' @export
#'
#' @examples
write_fasta <- function(x, file = "", compress = FALSE, quiet=FALSE) {
  if(stringr::str_detect(file, "\\.gz$") & !compress){
    compress <- TRUE
    if(!quiet) message(".gz detected in filename, compressing output file")
  }
  if (!is.null(dim(x))) {
    x <- as.list(as.data.frame(t(unclass(x))))
  }
  if (inherits(x, "DNAbin")) {
    tmp <- DNAbin2char(x)
  } else if (is.list(x)) {
    if (length(x[[1]] == 1)) {
      tmp <- unlist(x, use.names = TRUE)
    }
    else {
      tmp <- sapply(x, paste0, collapse = "")
    }
  } else {
    tmp <- x
  }
  reslen <- 2 * length(tmp)
  res <- character(reslen)
  res[seq(1, reslen, by = 2)] <- paste0(">", names(tmp))
  res[seq(2, reslen, by = 2)] <- tmp
  if(!file == ""){
    f <- if(compress){
        gzfile(file, "w")
      } else {
        file(file, "w")
    }
    writeLines(res, f)
    close(f)
  } else {
    writeLines(res)
  }
  if(!quiet){message("Wrote ", length(tmp), " seqeuences to ", file)}
  invisible(NULL)
}

#' Concatenate DNAbin objects while preserving attributes.
#' This function joins two or more \code{DNAbin} objects, retaining any
#'   attributes whose lengths match those of the input objects (e.g. "species",
#'   "lineage" and/or "taxID" attributes).
#' @param ... \code{DNAbin} objects, or list of DNAbins to be concatenated.
#'
#' @return
#' @export
#'
#' @examples
concat_DNAbin <- function(...){
  dots <- list(...)
  if(is.list(dots[[1]]) & length(dots) == 1){
    dots <- dots[[1]]
  }
  to_remove <- sapply(dots, is.null)
  dots <- dots[!to_remove]
  nlsts <- length(dots)
  DNA <- any(sapply(dots, class) == "DNAbin")
  if(nlsts == 0) return(NULL)
  if(nlsts == 1) return(dots[[1]])
  islist <- sapply(dots, is.list)
  for(i in which(!islist)){
    print(i)
    tmpattr <- attributes(dots[[i]])
    attributes(dots[[i]]) <- NULL
    dots[[i]] <- list(dots[[i]])

    print(tmpattr)
    attributes(dots[[i]]) <- tmpattr
  }
  names(dots) <- NULL #Remove name of list object to prevent name concatenation by unclass
  dots <- lapply(dots, unclass)
  findattr <- function(x){
    names(attributes(x))[sapply(attributes(x), length) == length(x)]
  }
  attrlist <- lapply(dots, findattr)
  ual <- unique(unlist(attrlist, use.names = FALSE))
  validattrs <- ual[sapply(ual, function(e) all(sapply(attrlist, function(g) e %in% g)))]
  validattrs <- validattrs[!validattrs %in% c("names", "class")]
  res <- unlist(dots, recursive = FALSE, use.names = TRUE)
  for(i in validattrs){
    attr(res, i) <- unlist(lapply(dots, attr, i), use.names = FALSE)
  }
  class(res) <- "DNAbin"
  return(res)
}



# .findExecutable ---------------------------------------------------------


#' Find executable
#'
#' @param exe the name of the executable to find
#' @param quiet Whether errors should be printed to console
#'
#' @return
#' @export
#'
.findExecutable <- function(exe, quiet=FALSE) {
  path <- Sys.which(exe)
  if(all(path=="")) {
    if(!quiet) stop("Executable for ", paste(exe, collapse=" or "), " not found! Please make sure that the software is correctly installed and, if necessary, path variables are set.", call.=FALSE)
    return(character(0))
  }

  path[which(path!="")[1]]
}


# Install BLAST -----------------------------------------------------------

#' Install BLAST
#'
#' @param url (Optional) Default will search for the latest version
#' URL to retrieve BLAST version from.
#' @param dest_dir (Optional)  Default "bin"
#' Directory to install BLAST within.
#' @param force Whether existing installs should be forcefully overwritten
#'
#' @return
#' @export
#' @import stringr
#' @importFrom RCurl getURL
#' @importFrom httr GET
#' @importFrom httr write_disk
#' @importFrom utils untar
#'
blast_install <- function(url, dest_dir = "bin", force = FALSE) {

  # get start time
  time <- Sys.time()
  # get OS
  localos <- Sys.info()["sysname"]

  if (missing(url)) {
    # find the latest version of BLAST
    url <- 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/'
    filenames <- RCurl::getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE) %>%
      stringr::str_split("\r*\n") %>%
      unlist()

    if(localos == "Windows"){
    url <- filenames[stringr::str_detect(filenames,"win64.tar.gz$")] %>%
      paste0(url,.)
    } else if(localos == "Darwin"){
      url <- filenames[stringr::str_detect(filenames,"macosx.tar.gz$")] %>%
        paste0(url,.)
    } else if(localos == "unix"){
      url <- filenames[stringr::str_detect(filenames,"linux.tar.gz$")] %>%
        paste0(url,.)
    }

  }

  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir) # Create first directory
  }

  blast_version <- basename(url) %>% stringr::str_replace("(-x64)(.*?)(?=$)", "")
  if (dir.exists(paste0(dest_dir, "/", blast_version)) && force == FALSE) {
    message("Skipped as BLAST already exists in directory, to overwrite set force to TRUE")
    return(NULL)
  } else  if (dir.exists(paste0(dest_dir, "/", blast_version)) && force == TRUE) {
    unlink(paste0(dest_dir, "/", blast_version), recursive = TRUE) # Remove old version
  }

  destfile <- file.path(dest_dir, basename(url))
  if (exists(destfile)) {
    file.remove(destfile) # Remove old zip file
  }

  #Download and unzip
  httr::GET(url, httr::write_disk(destfile, overwrite=TRUE))
  utils::untar(destfile, exdir = dest_dir)
  file.remove(destfile)

  #Set new $Paths variable for mac & linux
  if(localos == "Darwin" | localos == "unix"){
    old_path <- Sys.getenv("PATH")
    install_path <- list.dirs(dest_dir, full.names = TRUE)[str_detect(list.dirs(dest_dir, full.names = TRUE),"/bin$")]
    Sys.setenv(PATH = paste(old_path, normalizePath(install_path), sep = ":"))
  }

  time <- Sys.time() - time
  message(paste0("Downloaded ", blast_version, " in ", format(time, digits = 2)))
}


# Make Blast DB -----------------------------------------------------------

#' Make blast Database
#'
#' @param file (Required) A fasta file to create a database from.
#' @param dbtype (Optional) Molecule type of database, accepts "nucl" for nucleotide or "prot" for protein.
#' @param args (Optional) Extra arguments passed to BLAST
#' @param quiet (Optional) Whether progress should be printed to console, default is FALSE
#' @param remove_gaps (Optional) Whether gaps should be removed from the fasta file. Note that makeblastdb can fail if there are too many gaps in the sequence.
#'
#' @return
#' @export
#' @import stringr
#' @import ape
#' @importFrom R.utils gunzip
#'
make_blast_db <- function (file, dbtype = "nucl", args = NULL, quiet = FALSE, remove_gaps=TRUE) {
  time <- Sys.time() # get time
  .findExecutable("makeblastdb") # Check blast is installed
  if (is.null(args)){args <- ""}

  # Create a temp file for new blast DB
  tmpfile <- tempfile()

  if(remove_gaps) {
    seqs <- ape::del.gaps(ape::read.FASTA(file))
    write_fasta(seqs, tmpfile, compress = FALSE)
  } else {
    if (stringr::str_detect(file, ".gz")) {
      message("Unzipping file")
      R.utils::gunzip(filename=file, destname=tmpfile, remove=FALSE, overwrite=TRUE)
    } else {
      file.copy(file, tmpfile)
    }
  }

  # Run makeblastdb executable
  results <- system2(command = .findExecutable("makeblastdb"),
                     args = c("-in", tmpfile, "-dbtype", dbtype, args),
                     wait = TRUE,
                     stdout = TRUE)
  time <- Sys.time() - time
  if (!quiet) (message(paste0("made BLAST DB in ", format(time, digits = 2))))
  return(tmpfile)
}


#' Show BLAST parameters
#'
#' @param type (Required) Which BLAST function to display help page for
#'
#' @return
#' @export
#'
blast_params <- function(type = "blastn") {
  system(paste(.findExecutable(c(type)), "-help"))
}


# BLAST -------------------------------------------------------------------

#' Run BLAST search
#'
#' @param query (Required) Query sequence. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' @param db (Required) Reference sequences to conduct search against. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' If DNAbin, DNAStringSet or character string is provided, a temporary fasta file is used to construct BLAST database.
#' If db is set to "remote", this will conduct a search against NCBI nucleotide database.
#' @param type (Required) type of search to conduct, default 'blastn'
#' @param evalue (Required) Minimum evalue from search
#' @param output_format The output format to be returned.
#'  Default is tabular, which returns 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcov
#' See https://www.ncbi.nlm.nih.gov/books/NBK279684/ for custom formatting information
#' @param args (Optional) Extra arguments passed to BLAST
#' @param ungapped Whether ungapped alignment should be conducted. Default is FALSE.
#' @param quiet (Optional) Whether progress should be printed to console, default is FALSE
#' @param multithread Whether multithreading should be used, if TRUE the number of cores will be automatically detected, or provided a numeric vector to manually set the number of cores to use
#' @param remove_db_gaps Whether gaps should be removed from the fasta file used for the database. Note that makeblastdb can fail if there are too many gaps in the sequence.
#'
#' @return
#' @export
#' @import stringr
#' @import dplyr
#' @import future
#' @importFrom tibble enframe
#' @importFrom Biostrings writeXStringSet
#' @importFrom Biostrings DNA_ALPHABET
#' @importFrom ape as.DNAbin
#' @importFrom ape del.gaps
#' @importFrom stringr str_remove
#' @importFrom stringr str_extract
#' @importFrom stringr str_split_fixed
#' @importFrom stringr str_to_upper
#' @importFrom stringr str_split
#' @importFrom future availableCores
blast <- function (query, db, type="blastn", evalue = 1e-6,
                   output_format = "tabular", args=NULL, ungapped=FALSE,
                   quiet=FALSE, multithread=FALSE, remove_db_gaps = TRUE){
  time <- Sys.time() # get time
  # Create temp files
  tmp <- tempdir()
  tmpquery <- paste0(tmp, "/tmpquery.fa")
  tmpdb <- paste0(tmp, "/tmpquery.fa")

  .findExecutable(type) # check blast is installed

  #Setup multithreading
  ncores <- future::availableCores() -1
  if((isTRUE(multithread) | is.numeric(multithread) & multithread > 1) & db=="remote"){
    stop("Multithreading must be set to false for remote searches")
  } else if(isTRUE(multithread)){
    cores <- ncores
    if(!quiet){message("Multithreading with ", cores, " cores")}
  } else if (is.numeric(multithread) & multithread > 1){
    cores <- multithread
    if(cores > ncores){
      cores <- ncores
      message("Warning: the value provided to multithread is higher than the number of cores, using ", cores, " cores instead")
    }
    if(!quiet){message("Multithreading with ", cores, " cores")}
  } else if(isFALSE(multithread) | multithread==1){
    cores <- 1
  } else (
    stop("Multithread must be a logical or numeric vector of the numbers of cores to use")
  )
  nthreads <- ifelse(cores > 1, paste0("-num_threads ", unname(cores)), "")

  # Check outfmt
  if(output_format=="tabular"){
    #Custom tabular format
    parsecols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                   "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "qcovhsp")
    outfmt <- paste0("\"",6," ", paste(parsecols, collapse=" "),"\"")
  } else if(is.numeric(output_format)){
    outfmt <- output_format
    parsecols <- NULL
  } else if (!is.na(stringr::str_extract(output_format, "."))){
    outfmt <- paste0("\"",output_format,"\"")
    parsecols <- output_format %>%
      stringr::str_remove("^..") %>%
      stringr::str_split_fixed("\ ", n=Inf) %>%
      as.character()
  }

  # Database
  if(db=="remote"){
    db <- "nt"
    remote <- "-remote"
  } else if(inherits(db, "DNAbin")){
    if (!quiet) { message("Database input is DNAbin: Creating temporary blast database") }
    write_fasta(db, tmpdb)
    db <- make_blast_db(tmpdb, remove_gaps = remove_db_gaps)
    remote <- ""
  } else if (inherits(db, "DNAString") | inherits(db, "DNAStringSet")){
    if (!quiet) { message("Database input is DNAStringSet: Creating temporary blast database") }
    Biostrings::writeXStringSet(db, tmpdb)
    db <- make_blast_db(tmpdb, remove_gaps = remove_db_gaps)
    remote <- ""
  } else if (inherits(db, "character") &&  all(stringr::str_to_upper(stringr::str_split(db,"")[[1]]) %in% Biostrings::DNA_ALPHABET)) { # Handle text input
    if (!quiet) { message("Database input is character string: Creating temporary blast database") }
    if (nchar(db[1]) == 1) {db <- paste0(db, collapse = "")}
    db <- char2DNAbin(db)
    write_fasta(db, tmpdb)
    db <- make_blast_db(tmpdb, remove_gaps = remove_db_gaps)
    remote <- ""
  } else if (inherits(db, "character") &&  file.exists(file.path(db))){ # Handle filename
    db <- make_blast_db(db, remove_gaps = remove_db_gaps)
    remote <- ""
  } else {
    stop("Invalid BLAST database")
  }

  # Query
  if(inherits(query, "DNAbin")){
    if (!quiet) { message("Query is DNAbin: Creating temporary fasta file") }
    query <- ape::del.gaps(query)
    write_fasta(query, tmpquery)
    input <- tmpquery
  } else if (inherits(query, "DNAString") | inherits(query, "DNAStringSet")){
    if (!quiet) { message("Query is DNAString: Creating temporary fasta file") }
    query <- ape::as.DNAbin(query)
    query <- ape::del.gaps(query)
    write_fasta(query, tmpquery)
    input <- tmpquery
  }else if (inherits(query, "character") &&  all(stringr::str_to_upper(stringr::str_split(query,"")[[1]]) %in% Biostrings::DNA_ALPHABET)) { # Handle text query
    if (!quiet) { message("Query is character string: Creating temporary fasta file") }
    if (nchar(query[1]) == 1) {query <- paste0(query, collapse = "")}
    query <- char2DNAbin(query)
    query <- ape::del.gaps(query)
    write_fasta(query, tmpquery)
    input <- tmpquery
  } else if (inherits(query, "character") &&  file.exists(file.path(query))){ # Handle filenames
    input <- query
  }else {
    stop("Invalid BLAST query")
  }

  #Setup ungapped
  ungapped <- ifelse(ungapped, "-ungapped ","")

  #  Conduct BLAST search
  if (!quiet) { message("Running BLASTn query: ", paste(c("-db", db,
                                                          "-query", input,
                                                          remote,
                                                          "-outfmt", outfmt,
                                                          "-evalue", evalue,
                                                          args,
                                                          nthreads,
                                                          ungapped), collapse=" "))}
  results <- system2(command = .findExecutable(type),
                     args = c("-db", db,
                              "-query", input,
                              remote,
                              "-outfmt", outfmt,
                              "-evalue", evalue,
                              args,
                              nthreads,
                              ungapped),
                     wait = TRUE,
                     stdout = TRUE)

  # Check for errors and stop if detected
  if(any(stringr::str_detect(results, "Error:"))){
    err_to_print <- results[stringr::str_detect(results, "Error:")]
    stop(err_to_print[1])
  }

  # Remove the error message about nucleotides in title
  results <- results[!stringr::str_detect(results, "Title ends with at least 20 valid nucleotide characters")]

  # Parse BLAST results
  if(!is.null(parsecols)){
    out <- results %>%
      tibble::enframe() %>%
      tidyr::separate(col = value, into = parsecols,  sep = "\t", convert = TRUE)
  } else{
    message("Warning, result parsing is currently only supported for output_format = 'tabular', returning raw results")
    out <- results %>%
      tibble::enframe()
  }
  time <- Sys.time() - time
  if (!quiet) (message(paste0("finished BLAST in ", format(time, digits = 2))))

  # Clean up files
  if(file.exists(tmpdb)){file.remove(tmpdb)}
  if(file.exists(tmpquery)){file.remove(tmpquery)}
  return(out)
}


# BLAST_top_hit -----------------------------------------------------------

#' BLAST Top Hit
#'
#' @description Conduct BLAST search and return top hit
#' @param query (Required) Query sequence. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' @param db (Required) Reference sequences to conduct search against. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' If DNAbin, DNAStringSet or character string is provided, a temporary fasta file is used to construct BLAST database
#' @param type (Required) type of search to conduct, default 'blastn'
#' @param identity (Required) Minimum percent identity cutoff. Note that this is calculated using all alignments for each query-subject match.
#' @param coverage (Required) Minimum percent query coverage cutoff. Note that this is calculated using all alignments for each query-subject match.
#' @param evalue (Required) Minimum expect value (E) for saving hits
#' @param max_target_seqs (Required) Number of aligned sequences to keep. Even if you are only looking for 1 top hit keep this higher for calculations to perform properly.
#' @param max_hsp (Required) Maximum number of HSPs (alignments) to keep for any single query-subject pair.
#' @param ranks (Required) The taxonomic ranks contained in the fasta headers
#' @param delim (Required) The delimiter between taxonomic ranks in fasta headers
#' @param tie How to handle ties in top hit results. Options are to break ties by selecting the first hit (Default), or return all tied hits.
#' @param args (Optional) Extra arguments passed to BLAST
#' @param quiet (Optional) Whether progress should be printed to console, default is FALSE
#' @param remove_db_gaps Whether gaps should be removed from the fasta file used for the database. Note that makeblastdb can fail if there are too many gaps in the sequence.
#'
#' @return
#' @export
#' @import dplyr
#' @import IRanges
#' @importFrom tidyr separate
#'
blast_top_hit <- function(query, db, type="blastn",
                          identity=95, coverage=95, evalue=1e06, max_target_seqs=5, max_hsp=5,
                          ranks=c("Kingdom", "Phylum","Class", "Order", "Family", "Genus", "Species"), delim=";",
                          resolve_ties="first", args=NULL, remove_db_gaps = TRUE, quiet=FALSE,...){

  # set input filters in advance to speed up blast
  args <- paste("-perc_identity", identity, "-max_target_seqs", max_target_seqs, "-max_hsps", max_hsp, args)


  #Conduct BLAST
  result <- blast(query=query, type=type, db=db,
                  evalue = evalue,
                  args=args,
                  output_format = '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovs',
                  remove_db_gaps = remove_db_gaps) %>%
    dplyr::filter(!is.na(sseqid))

  # Check if ranks set up correctly
  if(!length(unlist(str_split(result$stitle[1], ";"))) == length(c("acc", ranks))){
    stop("number of ranks does not match database!")
  }

  #Subset to top hit
  # Blast may only align a portion of the sequence leaving gaps at the end
  # The full_pident calculation aims to correct this, calculating the percentage identity across full length of query sequence
  # Full-length pid modified from: https://github.com/McMahonLab/TaxAss/blob/master/tax-scripts/calc_full_length_pident.R
  top_hit <- result %>%
    dplyr::mutate(q_align = qend - qstart + 1, # Length of alignment between query and reference
                  q_len_adj = ifelse(qlen > slen, slen, qlen) # Handle case when query is longer than subject
    ) %>%
    dplyr::mutate(full_pident = (pident * length)/(length - q_align + q_len_adj)) %>%
    dplyr::group_by(qseqid, sseqid, stitle, qstart, qend, length, full_pident)  %>%
    dplyr::slice(1)%>% # Handle rare cases where there are multiple identical scoring matches within a single subject sequence (I.e multiple 16s copies)
    dplyr::ungroup() %>%
    dplyr::group_by(qseqid, sseqid, stitle) %>%
    dplyr::group_modify(~{
      # Handle cases where there are multiple matches overlapping the same segment of a subject by picking the highest scoring hit
      # This is common when paired end reads that do not overlap are joined together (i.e. with concat_unmerged in freyr)
      if(nrow(.x) > 1){ # Dont check cases with single unique hits

        # Check if hits overlap the same region
        # setup the IRanges object from the input qstart and qend
        ir <- IRanges::IRanges(as.numeric(.x$qstart), as.numeric(.x$qend), names = .x$name)

        # find which hit ids overlap with each other
        ovrlp <- IRanges::findOverlaps(ir, drop.self = TRUE, drop.redundant = TRUE)

        # store id indices for further use
        hit1 <- queryHits(ovrlp)
        hit2 <- subjectHits(ovrlp)

        # width of overlaps between ids
        widths <- width(pintersect(ir[hit1], ir[hit2])) - 1

        # result
        overlaps <- data.frame(id1 = names(ir)[hit1], id2 = names(ir)[hit2], widths)

        # if the multiple hits are overlapping, get the best hit - otherwise leave them as they will have been handled correctly when summing full_pident
        if(nrow(overlaps) > 0){
          newdf <- list()
          for (i in 1:nrow(overlaps)){
            newdf[[i]] <- .x %>%
              filter(as.character(name) %in% c(overlaps$id1[i], overlaps$id2[i]))%>%
              dplyr::top_n(1, bitscore) %>%
              dplyr::top_n(1, pident) %>%
              dplyr::top_n(1, qcovs)
          }
          return(newdf %>%
                   bind_rows() %>%
                   distinct())
        } else {
          return(.x)
        }
      } else {
        return(.x)
      }
    })%>%
    dplyr::summarise(pident = sum(full_pident), # Combine hit stats for multiple discontiguous matches
                     qcovs = unique(qcovs),
                     max_score = max(bitscore),
                     total_score = sum(bitscore),
                     evalue = min(evalue)) %>%
    dplyr::filter(pident > identity, qcovs > coverage) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(qseqid) %>%
    dplyr::top_n(1, total_score) %>%
    dplyr::top_n(1, max_score) %>%
    dplyr::top_n(1, qcovs) %>%
    dplyr::top_n(1, pident) %>%
    tidyr::separate(stitle, c("acc", ranks), delim, remove = TRUE)

  if(resolve_ties == "first"){
    top_hit <- top_hit %>%
      dplyr::mutate(row_n = dplyr::row_number()) %>%
      dplyr::top_n(1, row_n) %>% # Break ties by position
      dplyr::select(-row_n) %>%
      dplyr::ungroup()
  } else if(resolve_ties == "all"){ # Return all ties
    top_hit <- top_hit %>%
      dplyr::ungroup()
  }
  return(top_hit)
}


# BLAST assign species ----------------------------------------------------


#' Assign species using BLAST
#'
#' @description This is to be used alongside a hierarchial classifier such as IDTAXA or RDP to assign additional species level matches.
#' This is designed to be a more flexible version of dada2's assignSpecies function
#' @param query (Required) Query sequence. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' @param db (Required) Reference sequences to conduct search against. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' If DNAbin, DNAStringSet or character string is provided, a temporary fasta file is used to construct BLAST database
#' @param type (Required) type of search to conduct, default 'blastn'
#' @param identity (Required) Minimum percent identity cutoff. Note that this is calculated using all alignments for each query-subject match.
#' @param coverage (Required) Minimum percent query coverage cutoff. Note that this is calculated using all alignments for each query-subject match.
#' @param evalue (Required) Minimum expect value (E) for saving hits
#' @param max_target_seqs (Required) Number of aligned sequences to keep. Even if you are only looking for 1 top hit keep this higher for calculations to perform properly.
#' @param max_hsp (Required) Maximum number of HSPs (alignments) to keep for any single query-subject pair.
#' @param ranks (Required) The taxonomic ranks contained in the fasta headers
#' @param delim (Required) The delimiter between taxonomic ranks in fasta headers
#' @param args (Optional) Extra arguments passed to BLAST
#' @param quiet (Optional) Whether progress should be printed to console, default is FALSE
#' @param remove_db_gaps Whether gaps should be removed from the fasta file used for the database. Note that makeblastdb can fail if there are too many gaps in the sequence.
#'
#' @return
#' @export
#' @import dplyr
#' @import stringr
#' @importFrom tidyr separate
#'
blast_assign_species <- function(query, db, type="blastn",
                                 identity=97, coverage=95, evalue=1e06, max_target_seqs=5, max_hsp=5,
                                 ranks=c("Kingdom", "Phylum","Class", "Order", "Family", "Genus", "Species"), delim=";",
                                 args=NULL, quiet=FALSE, remove_db_gaps=TRUE){

  #Check input contains species and genus
  if(!any(tolower(ranks) %in% c("species", "genus"))){
    stop("Ranks must include Genus and Species")
  }

  #Conduct BLAST
  top_hit <- blast_top_hit(query = query, db = db, type=type,
                           identity=identity, coverage=coverage, evalue=evalue, max_target_seqs=max_target_seqs, max_hsp=max_hsp,
                           ranks=ranks, delim=delim, resolve_ties="all", args=args, quiet=quiet, remove_db_gaps=remove_db_gaps) %>%
    dplyr::filter(!is.na(Species))

  out <- top_hit %>%
    dplyr::group_by(qseqid) %>%
    dplyr::mutate(spp = Species %>% stringr::str_remove("^.* ")) %>%
    dplyr::reframe(spp = paste(sort(unique(spp)), collapse = "/"), Genus, pident, qcovs, max_score, total_score, evalue) %>%
    dplyr::mutate(binomial = paste(Genus, spp)) %>%
    dplyr::distinct() %>%
    dplyr::group_by(qseqid) %>%
    dplyr::add_tally() %>%
    dplyr::mutate(binomial =  dplyr::case_when( #Leave unassigned if conflicted at genus level
      n > 1 ~ as.character(NA),
      n == 1 ~ binomial
    )) %>%
    dplyr::select(OTU = qseqid, Genus, Species = binomial, pident, qcovs, max_score, total_score, evalue)

  return(out)
}

# Genbank functions -----------------------------------------------
# Some useful entrez queries
#' all `[filter]` 	Retrieves everthing
#' Specified `[property]` 	Formal binomial and trinomial
#' at or below species level `[property]`
#' family `[rank]` 	Rank-based query
#' taxonomy genome `[filter]` 	Taxa with a direct link to a genome sequence
#' 2009/10/21:2020 `[date]` 	Date-bounded query
#' mammalia `[subtree]` 	All taxa within the Mammalia
#' extinct `[property]` 	Extinct organisms
#' Terminal `[property]` 	Terminal nodes in the tree
#' loprovencyclife `[filter]` 	Entries with LinkOut links to the Encyclopedia of Life
#'

#' Fetch sequences from genbank
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param database The database to download from. For NCBI GenBank this currently onlt accepts the arguments 'nuccore' or 'genbank' which is an alias for nuccore.
#' @param marker The barcode marker used as a search term for the database.
#' If this is set to "mitochondria" or "mitochondrion" it will download full mitochondrial genomes. If set to "genome" it will download entire genomes only.
#' @param output The output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (Accession;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (Accession;Genus_species),
#' "bold" for BOLD taxonomic ID only (Accession;BoldTaxID),
#' "gb" for genbank taxonomic ID (Accession|GBTaxID),
#' "gb-binom" which outputs genbank taxonomic ID's and Genus species binomials, translating BOLD taxonomic ID's to genbank in the process (Accession|GBTaxID;Genus_species)
#' or "standard" which outputs the default format for each database. For bold this is `sampleid|species name|markercode|genbankid`
#' @param min_length The minimum length of sequences to download
#' @param max_length The maximum length of sequences to download
#' @param subsample (Numeric) return a random subsample of sequences from the search.
#' @param chunk_size Split up the query into chunks of this size to avoid overloading API servers. if left NULL, the default will be 300
#' @param db a database file generated using ` get_ncbi_taxonomy()`. Generated automatically if NULL.
#' @param multithread Whether multithreading should be used, if TRUE the number of cores will be automatically detected, or provided a numeric vector to manually set the number of cores to use
#' @param quiet Whether progress should be printed to the console.
#' @param progress A logical, for whether or not to print a progress bar when multithread is true. Note, this will slow down processing.
#' @param retry_attempt The number of query attempts in case of query failure due to poor internet connection.
#' @param retry_wait How long to wait between query attempts.
#'
#' @return
#'
#' @examples
fetch_genbank <- function(x, database = "nuccore", marker = c("COI[GENE]", "CO1[GENE]", "COX1[GENE]"), output = "h",
                           min_length = 1, max_length = 2000, subsample=FALSE, chunk_size=100, db=NULL,
                           multithread = FALSE, quiet = FALSE, progress=FALSE, retry_attempt=3, retry_wait=5) {
  if(!database %in% c("nuccore", "genbank")){
    stop("database is invalid: only nuccore and genbank is currently supported")
  }
  if (output == "bold"){
    stop("bold output is only valid for searching the bold database")
  }
  if (!output %in% c("standard", "h", "binom", "gb", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'standard', 'h','binom','gb' or 'gb-binom', see help page for more details"))
  }
  if (!any(tolower(marker) %in% c("mitochondria", "mitochondrion", "genome"))) {
    marker <- paste(marker, collapse=" OR ")
  }

  # Get NCBI taxonomy database if NCBI format outputs are desired
  if (is.null(db) && output %in% c("gb", "gb-binom", "h")) {
    db <- get_ncbi_taxonomy(include_synonyms = TRUE, force=FALSE)
  }

  if(database=="genbank"){
    database="nuccore"
  }
  # Setup multithreading
  setup_multithread(multithread, quiet=quiet)

  # Main function
  if (!tolower(marker) %in% c("mitochondria", "mitochondrion", "genome")) {
    searchQ <- paste("(", x, "[ORGN])", " AND (", paste(c(marker), collapse = " OR "), ") AND ", min_length, ":", max_length, "[Sequence Length]", sep = "")
  } else if (tolower(marker) %in% c("mitochondria", "mitochondrion")) {
    message(paste0("Input marker is ", marker, ", Downloading full mitochondrial genomes"))
    searchQ <- paste("(", x, "[ORGN])", " AND mitochondrion[filter] AND genome", sep = "")
  } else if (tolower(marker) %in% c("genome")){
    message(paste0("Input marker is ", marker, ", Downloading full genomes"))
    searchQ <- paste("(", x, "[ORGN])", " AND complete genome[title]", sep = "")
  }
  if(!quiet){message("Searching genbank with query:", searchQ)}

  search_results <- rentrez::entrez_search(db = database, term = searchQ, retmax = 9999999, use_history = TRUE)
  if (search_results$count > 0 & !is.na(search_results$count)) {

    if (is.numeric(subsample)){
      if (!quiet) {message(paste0(subsample, " of: ", search_results$count, " sequences to be downloaded for: ", searchQ))}
      ids <- sample(search_results$ids, subsample )
    } else{
      if (!quiet) {message(paste0(search_results$count, " Sequences to be downloaded for: ", searchQ))}
      ids <- search_results$ids
    }

    # Split query into chunks
    n <- length(ids)
    r <- rep(1:ceiling(n/chunk_size), each=chunk_size)[1:n]
    id_list <- split(ids, r) %>% magrittr::set_names(NULL)

    # Retrieve each chunk
    seqs <- furrr::future_map(
      id_list, read_genbank_chunk, quiet = FALSE, retry_attempt = retry_attempt, retry_wait = retry_wait, .progress = progress)

    failed <- id_list[!sapply(seqs, class)=="DNAbin"]
    seqs <- seqs[sapply(seqs, class)=="DNAbin"]
    seqs <- concat_DNAbin(seqs)

    # Parse attributes
    seq_spp <- attr(seqs, "species") %>% stringr::str_replace_all("_", " ")
    seq_acc <- attr(seqs, "acc") %>% stringr::str_remove(" .*$")

    if (output == "standard") { # Standard output
      names <- paste0(seq_acc, ";", names(seqs))
    } else if (output == "binom") {
      names <- paste0(seq_acc, ";", seq_spp)
    } else if (output == "h") { # Hierarchical output
      names <- tibble::enframe(seq_spp, name=NULL, value="tax_name") %>%
        dplyr::mutate(sampleid = seq_acc) %>%
        dplyr::left_join(db, by="tax_name") %>%
        tidyr::unite("names", c("sampleid", "superkingdom", "phylum", "class", "order", "family", "genus", "tax_name"), sep = ";")%>%
        dplyr::mutate(names = names %>%
                        stringr::str_replace_all(pattern = ";;", replacement = ";")) %>%
        dplyr::pull(names)
    } else if (output == "gb") { # Genbank taxID output
      names <- tibble::enframe(seq_spp, name=NULL, value="tax_name") %>%
        dplyr::mutate(sampleid = seq_acc) %>%
        dplyr::left_join(db, by="tax_name") %>%
        tidyr::unite("name", c("sampleid", "tax_id"), sep = "|") %>%
        dplyr::pull(name)
    } else if (output == "gb-binom") {
      names <- tibble::enframe(seq_spp, name=NULL, value="tax_name") %>%
        dplyr::mutate(sampleid = seq_acc) %>%
        dplyr::left_join(db, by="tax_name") %>%
        tidyr::unite("name", c("sampleid", "tax_id"), sep = "|") %>%
        tidyr::unite("name", c("name", "tax_name"), sep = ";") %>%
        dplyr::pull(name)
    }
    names(seqs) <- names %>% stringr::str_replace_all(" ", "_")
    out <- seqs
    attr(out, "failed") <- failed
  } else {
    if (!quiet)(message("No sequences available for query: ", searchQ))
    out <- NULL
  }
  # Explicitly close multisession workers
  future::plan(future::sequential)
  return(out)
}

#' Read genbank chunk
#'
#' @param gid a vector of GenBank ID's
#' @param quiet Whether progress should be printed to console.
#' @param retry_attempt The number of query attempts in case of query failure due to poor internet connection.
#' @param retry_wait How long to wait between query attempts
#'
#' @return
#'
#' @examples
read_genbank_chunk <- function(gid, quiet = FALSE, retry_attempt=3, retry_wait=5) {
  n_seqs <- length(gid)
  URL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", paste(gid, collapse = ","), "&rettype=gb&retmode=text", sep = "")
  gb <- NULL
  seqs <- NULL
  attempt <- 1
  # Download with error handling
  while((is.null(gb) | (length(seqs) < length(gid))) && attempt <= (retry_attempt+1)) {
    gb <- tryCatch({
      scan(file = URL, what = "", sep = "\n", quiet = TRUE)
    }, error = function(e){
      if (!quiet) {cat(paste("Failed attempt ", attempt,"\n"))}
      Sys.sleep(retry_wait)
      NULL
    })
    if(!is.null(gb)){
      seqs <- parse_gb(gb)
    }
    attempt <- attempt + 1
  }
  if(!length(seqs) == length(gid)){
    if(!quiet)warning("length of returned sequences does not match length of query")
  }
  out <- seqs
  if(!is.null(out)){
    attr(out, "query") <- gid
  }
  return(out)
}

#' Parse genbank flat files
#'
#' @param gb A genbank flat file
#'
#' @return
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_detect
#' @examples
parse_gb <- function(gb){
  # Truncate record at last // to avoid broken records
  gb <- tryCatch({
    gb[1:max(which(grepl("^//", gb)))]
  }, error = function(e){
    warning("Failed parsing")
    print(max(which(grepl("^//", gb))))
    return(NULL)
  })
  if(sum(grepl("ACCESSION", gb)) < 1){
    return(NULL)
  }
  start <- which(grepl("^ORIGIN", gb))
  stop <- which(grepl("^//", gb))

  # Check for malformed start and stops
  good_records <- stringr::str_detect(gb[stop-1], "[0-9] [a-z-]")
  stop <- stop[good_records]

  if (length(start) == length(stop)){
    n_seqs <- length(start)
  } else{
    writeLines(gb, "failed.gb")
    stop("Incorrect length of start and stop, dumped failing records to failed.gb")
  }
  seqs <- vector("character", length = n_seqs)
  for (l in 1:n_seqs){
    seqs[[l]] <- toupper(paste(stringr::str_remove_all(gb[(start[l]+1):(stop[l]-1)], "[^A-Za-z]"), collapse=""))
  }
  # Extract relevant metadata
  seq_acc <- gsub("+ACCESSION +", "", grep("ACCESSION", gb, value = TRUE))
  seq_defs <- gsub("+DEFINITION +", "", grep("DEFINITION", gb, value = TRUE))
  seq_spp <- gsub(" ", "_", gsub(" +ORGANISM +", "", grep(" +ORGANISM +", gb, value = TRUE)))

  names(seqs) <- seq_defs[good_records]
  seqs <- char2DNAbin(seqs)
  attr(seqs, "acc") <- seq_acc[good_records]
  attr(seqs, "species") <- seq_spp[good_records]
  return(seqs)
}

# BOLD functions --------------------------------------------------------

#' Fetch sequences from BOLD
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param marker the barcode marker used as a search term for the database
#' @param quiet Whether progress should be printed to the console.
#' @param output the output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (Accession;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (Accession;Genus_species),
#' "bold" for BOLD taxonomic ID only (Accession;BoldTaxID),
#' "gb" for genbank taxonomic ID (Accession|GBTaxID),
#' "gb-binom" which outputs genbank taxonomic ID's and Genus species binomials, translating BOLD taxonomic ID's to genbank in the process (Accession|GBTaxID;Genus_species)
#' or "standard" which outputs the default format for each database. For genbank this is the description field, and for bold this is `sampleid|species name|markercode|genbankid`
#' @param db (Optional) a database file generated using ` get_ncbi_taxonomy()` or ` get_ott_taxonomy()`
#' @param retry_attempt The number of query attempts in case of query failure due to poor internet connection.
#' @param retry_wait How long to wait between query attempts
#'
#' @import dplyr
#' @import stringr
#' @importFrom bold bold_seqspec
#' @importFrom tidyr unite
#' @importFrom tidyr drop_na
#' @importFrom methods is
#' @return
#'
#' @examples
fetch_bold <- function(x, marker = "COI-5P", output = "gb-binom", quiet = FALSE, db=NULL, retry_attempt=3, retry_wait=5) {
  # function setup
  if (!output %in% c("h", "standard", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'standard', 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }

  # Get NCBI taxonomy database if NCBI format outputs are desired
  if (is.null(db) && output %in% c("gb", "gb-binom")) {
    db <- get_ncbi_taxonomy(include_synonyms = TRUE, force=FALSE)
  }

  # Bold search
  data <- tryCatch({
    read_bold_chunk(taxon = x, quiet = FALSE, retry_attempt = retry_attempt, retry_wait = retry_wait)
  }, error = function(e){
    writeLines(x, "failed.txt")
    stop("An error during data download, wrote failed taxa to failed.txt")
  })

  bckup_seqs <- data
  if (length(data) >0 & methods::is(data, "data.frame")) {

    data <- tryCatch({
      data %>%
        #dplyr::na_if("") %>%
        dplyr::mutate(domain_name = "Eukaryota") %>%
        dplyr::filter(markercode == marker) %>% # Remove all sequences for unwanted markers
        dplyr::filter(!is.na(species_name), !is.na(markercode), !is.na(nucleotides))
    }, error = function(e){
      saveRDS(bckup_seqs, "bold_data.rds")
      stop("An error during first stage of data parsing, dumped intermediate files to bold_data.rds")
    })

    if (nrow(data) > 0) {

      data <- tryCatch({
        if (output == "standard") {
          data %>%
            dplyr::select(sampleid, species_name, markercode, genbank_accession, nucleotides) %>%
            dplyr::mutate(species_name = species_name %>%
                            stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                            trimws(which = "both")) %>%
            tidyr::unite("name", c("sampleid", "species_name", "markercode", "genbank_accession"), sep = "|")
        } else if (output == "h") {
          # Hierarchial output
          data %>%
            dplyr::select(sampleid, domain_name, phylum_name, class_name,
                          order_name, family_name, genus_name, species_name,nucleotides
            ) %>%
            tidyr::drop_na() %>%
            tidyr::unite("name", c(
              "sampleid", "domain_name",
              "phylum_name", "class_name",
              "order_name", "family_name",
              "genus_name", "species_name"
            ),
            sep = ";"
            ) %>%
            dplyr::mutate(name = name %>%
                            stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                            trimws(which = "both"))
        } else if (output == "binom") {
          # Binomial output
          data %>%
            dplyr::select(sampleid, species_name, nucleotides) %>%
            tidyr::drop_na() %>%
            tidyr::unite("name", c("sampleid", "species_name"), sep = ";") %>%
            dplyr::mutate(name = name %>%
                            stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                            trimws(which = "both"))
        } else if (output == "bold") {
          # BOLD taxID output
          data %>%
            dplyr::select(sampleid, species_taxID, nucleotides) %>%
            tidyr::drop_na() %>%
            tidyr::unite("name", c("sampleid", "species_taxID"), sep = "|") %>%
            dplyr::mutate(name = name %>%
                            stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                            trimws(which = "both"))

        } else if (output == "gb") {
          # Genbank taxID output
          data %>%
            dplyr::select(sampleid, species_name, nucleotides) %>%
            tidyr::drop_na() %>%
            dplyr::mutate(tax_name = trimws(species_name, which = "both")) %>%
            dplyr::left_join(db, by="tax_name")%>%
            dplyr::mutate(tax_name = tax_name %>%
                            stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                            trimws(which = "both")) %>%
            tidyr::unite("name", c("sampleid", "tax_id"), sep = "|")

        } else if (output == "gb-binom") {
          # genbank binomial output
          data %>%
            dplyr::select(sampleid, species_name, nucleotides) %>%
            tidyr::drop_na() %>%
            dplyr::mutate(tax_name = trimws(species_name, which = "both")) %>%
            dplyr::left_join(db, by="tax_name") %>%
            dplyr::mutate(tax_name = tax_name %>%
                            stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                            trimws(which = "both")) %>%
            tidyr::unite("name", c("sampleid", "tax_id"), sep = "|") %>%
            tidyr::unite("name", c("name", "tax_name"), sep = ";")
        }
      }, error = function(e){
        saveRDS(bckup_seqs, "bold_data.rds")
        stop("An error during second stage of data parsing, dumped intermediate files to bold_data.rds")
      })
      seqs <- char2DNAbin(data$nucleotides)
      names(seqs) <- data$name
      out <- seqs
    } else {
      if (!quiet)(message("No sequences available for query: ", x, " and marker: ", marker))
      out <- NULL
    }
  } else {
    if (!quiet)(message("No sequences available for query: ", x, " and marker: ", marker))
    out <- NULL
  }
  return(out)
}

#' read bold chunk
#'
#' @param taxon A taxon name to download sequences for
#' @param gid a vector of GenBank ID's
#' @param quiet Whether progress should be printed to console.
#' @param retry_attempt The number of query attempts in case of query failure due to poor internet connection.
#' @param retry_wait How long to wait between query attempts
#'
#' @return
#' @export
#'
#' @examples
read_bold_chunk <- function(taxon, quiet = FALSE, retry_attempt=3, retry_wait=5) {
  dl <- NULL
  out <- NULL
  attempt <- 1
  # Download with error handling
  while(is.null(dl) && attempt <= (retry_attempt+1)) {
    dl <- tryCatch({
      cli <- crul::HttpClient$new(url = 'https://v4.boldsystems.org/index.php/API_Public/combined')
      out <- cli$get(query = list(taxon = taxon, format = "tsv"))
      out$raise_for_status()
      if (grepl("html", out$response_headers$`content-type`)) {
        stop(out$parse("UTF-8"))
      }
      out
    }, error = function(e){
      if (!quiet) {cat(paste("Failed attempt ", attempt,"\n"))}
      Sys.sleep(retry_wait)
      NULL
    })
    if(!is.null(dl)){
      tt <- paste0(rawToChar(dl$content, multiple = TRUE), collapse = "")
      if (tt == "") {
        return(NULL)
      }
      Encoding(tt) <- "UTF-8"
      if (grepl("Fatal error", tt)) {
        # set dl to null again to loop through again
        dl <- NULL
      } else{
        out <- readr::read_tsv(tt)
      }
    }
    attempt <- attempt + 1
  }
  return(out)
}


#' Split bold query
#' This function recursively splits a bold taxonomic query until the amount of records returned is under chunk_size
#'
#' @param x The input taxonomic query
#' @param chunk_size The maximum amount of records to return per query
#' @param split_if_under whether to split the query down 1 level even if already below chunksize. This can be useful for making suitable queries for multithread searching.
#' @param quiet Whether progress should be printed to the console.
#'
#' @return
#' @export
#'
#' @import dplyr
#' @importFrom taxize downstream
#' @importFrom bold bold_stats
#' @importFrom methods as
#'
split_bold_query <- function(x, chunk_size=100000, split_if_under = FALSE, quiet=FALSE){

  rcds <- bold::bold_stats(x, dataType = "overview") %>%
    unlist()
  out <- character()
  if(rcds["total_records"] > chunk_size){
    do_split <- TRUE
    if(!quiet) {message("Found over ", chunk_size, " (chunk_size value) BOLD records for ", x, ", searching for lower taxonomic ranks to reduce query size")}
  } else if(split_if_under & (rcds["species.count"] > 1)) {
    do_split <- TRUE
  } else {
    do_split <- FALSE
  }
  if(isTRUE(do_split)){
    while(length(x) > 0){
      rcds <- purrr::map(x, ~{
        .x %>%
          bold::bold_stats(dataType = "overview") %>%
          unlist()
      }) %>%
        dplyr::bind_rows() %>%
        dplyr::select("order.count", "family.count", "genus.count", "species.count") %>%
        dplyr::select_if(colSums(.) > nrow(.))

      downstream2 <- stringr::str_remove(colnames(rcds[1]), ".count")
      lower_ranks <- taxize::downstream(x, db = "bold", downto = downstream2) %>%
        methods::as("list") %>%
        dplyr::bind_rows() %>%
        dplyr::filter(rank == stringr::str_to_lower(!!downstream2)) %>%
        dplyr::pull(name)

      # Check if any are still over (Only if rank is higher than species)
      if (!downstream2 == "species"){
        rcds2 <- purrr::map(lower_ranks, ~{
          .x %>%
            bold::bold_stats(dataType = "overview") %>%
            unlist()
        }) %>%
          purrr::set_names(lower_ranks) %>%
          dplyr::bind_rows(.id="name")
        # Add successfully resolved to out
        out <- c(out, rcds2 %>%
                   dplyr::filter(total_records < chunk_size)%>%
                   dplyr::pull(name))
        # Repeat on unresolved
        x <- rcds2 %>%
          dplyr::filter(total_records > chunk_size) %>%
          dplyr::pull(name)
      } else {
        out <- c(out, lower_ranks)
        x <- character()
      }
    }
  } else(
    out <- x
  )
  return(out)
}

# wrapper -----------------------------------------------------------------

#' fetch_seqs function
#'
#' @param x A taxon name or vector of taxon names to download sequences for.
#' @param database The database to download from. For NCBI GenBank this currently onlt accepts the arguments 'nuccore' or 'genbank' which is an alias for nuccore.
#' Alternatively sequences can be downloaded from the Barcode of Life Data System (BOLD) using 'bold'
#' @param marker The barcode marker used as a search term for the database. If you are targetting a gene, adding a suffix \[GENE\] will increase the search selectivity.
#' The default for Genbank is 'COI\[GENE\] OR COX1\[GENE\] OR COXI\[GENE\]', while the default for BOLD is 'COI-5P'.
#' If this is set to "mitochondria" and database is 'nuccore', or 'genbank'it will download mitochondrial genomes only.
#' If this is set to "genome" and database is 'nuccore', or 'genbank'it will download complete genome sequences only.
#' @param output The output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (Accession;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (Accession;Genus_species),
#' "bold" for BOLD taxonomic ID only (Accession;BoldTaxID),
#' "gb" for genbank taxonomic ID (Accession|GBTaxID),
#' "gb-binom" which outputs genbank taxonomic ID's and Genus species binomials, translating BOLD taxonomic ID's to genbank in the process (Accession|GBTaxID;Genus_species)
#' or "standard" which outputs the default format for each database. For bold this is `sampleid|species name|markercode|genbankid`
#' @param min_length The maximum length of the query sequence to return. Default 1.
#' @param max_length The maximum length of the query sequence to return.
#' This can be useful for ensuring no off-target sequences are returned. Default 2000.
#' @param subsample (Numeric) return a random subsample of sequences from the search.
#' @param chunk_size Split up the queries made (for genbank), or returned records(for BOLD) into chunks of this size to avoid overloading API servers.
#' if left NULL, the default for genbank searches will be 10,000 for regular queries, 1,000 if marker is "mitochondria", and 1 if marker is "genome"
#' For BOLD queries the default is 100,000 returned records
#' @param db a database file generated using ` get_ncbi_taxonomy()`. Generated automatically if NULL.
#' @param multithread Whether multithreading should be used, if TRUE the number of cores will be automatically detected, or provided a numeric vector to manually set the number of cores to use
#' @param quiet Whether progress should be printed to the console.
#' @param progress A logical, for whether or not to print a progress bar when multithread is true. Note, this will slow down processing.
#' @param retry_attempt The number of query attempts in case of query failure due to poor internet connection.
#' @param retry_wait How long to wait between query attempts.
#'
#' @import dplyr
#' @import stringr
#' @import purrr
#' @import future
#' @import furrr
#' @importFrom taxize downstream
#' @importFrom methods as
#'
#' @return
#' @export
#'
#' @examples
fetch_seqs <- function(x, database, marker = NULL, output = "gb-binom",
                        min_length = 1, max_length = 2000, subsample=FALSE,
                        chunk_size=NULL, db=NULL, multithread = FALSE,
                        quiet = FALSE, progress=FALSE, retry_attempt=3, retry_wait=5) {
  # function setup
  time <- Sys.time() # get time

  database <- tolower(database)
  if(!database %in% c("nuccore", "genbank", "bold")) {
    stop("database is invalid. See help page for more details")
  }

  if (!output %in% c("standard", "h", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'standard', 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }

  #stop if subsample and BOLD is true
  if (is.numeric(subsample ) & database=="bold"){ stop("Subsampling is currently not supported for BOLD")}

  # Get NCBI taxonomy database if NCBI format outputs are desired
  if (is.null(db) & output %in% c("gb", "gb-binom")) {
    db <- get_ncbi_taxonomy(include_synonyms = TRUE, force=FALSE)
  }

  # Make sure x is unique
  x <- unique(x)

  # Genbank
  if (database %in% c("genbank", "nuccore")) {
    if(!quiet) {message("Downloading from genbank")}
    if(is.null(chunk_size)){
      chunk_size <- 100
    }
    res <- purrr::map(x, fetch_genbank, database = database, marker = marker,
                      output = output, min_length = min_length, max_length = max_length,
                      subsample=subsample, chunk_size=chunk_size, quiet = quiet, multithread=multithread,
                      db=db, retry_attempt = retry_attempt, retry_wait = retry_wait, progress = progress)
    res <- res[!sapply(res, is.null)]
    res <- concat_DNAbin(res)
    res
  } else if (database == "bold") {

    setup_multithread(multithread = multithread, quiet=quiet)

    # Split any very large queries to avoid overloading BOLD query
    if(is.null(chunk_size)){
      chunk_size <- 100000
    }
    # Also If multithreading and only 1 taxon provided, split main query into lower level taxonomy to create something to loop across
    bold_taxon <- furrr::future_map(
      x, split_bold_query, chunk_size=chunk_size, quiet=quiet,
      split_if_under = ifelse((isTRUE(multithread) | (is.numeric(multithread) & multithread > 1)) && length(x) == 1,
                              TRUE, FALSE), .options = furrr::furrr_options(seed = TRUE)) %>%
      unlist()

    if(!quiet) {message("Downloading ", length(bold_taxon)," taxa from BOLD")}
    res <- furrr::future_map(
      bold_taxon, fetch_bold, marker = marker, db=db, output = output,
      quiet = quiet, .progress = progress, .options = furrr::furrr_options(seed = TRUE))
    res <- res[!sapply(res, is.null)]
    res <- concat_DNAbin(res)
    res
  }

  # Explicitly close multisession workers
  future::plan(future::sequential)
  time <- Sys.time() - time
  if (!quiet) (message(paste0("Downloaded ", length(res), " ", x, " Sequences from ",database, " in ", format(time, digits = 2))))
  # Return results summary
  return(res)
}



# Update ------------------------------------------------------------------

#These functions are deprecated and are to be removed

#update_genbank <- function(x, fasta, database = "nuccore", marker = c("COI[GENE]", "CO1[GENE]", "COX1[GENE]"), quiet = FALSE, output = "h", suffix="updates",
#                           min_length = 1, max_length = 2000, chunk_size=300, out_dir = NULL,
#                           compress = FALSE, force=FALSE, multithread = FALSE, progress=FALSE){
#
#  # function setup
#  time <- Sys.time() # get time
#
#  if(!database %in% c("nuccore", "genbank")){
#    stop("database is invalid: only nuccore and genbank is currently supported")
#  }
#
#  #Check if all inputs exists
#  if (!any(file.exists(fasta))) {
#    stop("Not all fasta files provided can be found, check the diferectory you have provided")
#  }
#
#  if (!output %in% c("standard", "h", "binom", "gb", "bold", "gb-binom")) {
#    stop(paste0(output, " has to be one of: 'standard', 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
#  }
#  if (!any(tolower(marker) %in% c("mitochondria", "mitochondrion", "genome"))) {
#    marker <- paste(marker, collapse=" OR ")
#  }
#  #Define directories
#  if (is.null(out_dir)) {
#    out_dir <- database
#    if (!quiet) (message(paste0("No input out_dir given, saving output file to: ", out_dir)))
#  }
#  if (!dir.exists(out_dir)) {
#    dir.create(out_dir)
#  }
#
#  # Create output file
#  name <- marker %>%
#    stringr::str_remove_all(pattern="\\[GENE]") %>%
#    stringr::str_remove_all(pattern="OR ") %>%
#    stringr::str_replace_all(pattern=" ", replacement ="_")
#
#  if (compress) {
#    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, "_", suffix, ".fa.gz")
#  } else if (!compress) {
#    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, "_", suffix, ".fa")
#  }
#
#  #Check if output exists
#  if (file.exists(out_file) && force==TRUE) {
#    file.remove(out_file)
#    cat("", file=out_file)
#  } else if (file.exists(out_file) && force==FALSE){
#    stop(paste0(out_file, " exists, set force = TRUE to overwrite"))
#  }
#
#  if(database=="genbank"){
#    database="nuccore"
#  }
#
#  # Get accessions from existing fastas
#  current <- acc_from_fasta(fasta)
#  if(!quiet){message(length(current), " unique accessions in fastas")}
#
#  # Genbank Search
#  if (!tolower(marker) %in%c("mitochondria", "mitochondrion", "genome")) {
#    searchQ <- paste("(", x, "[ORGN])", " AND (", paste(c(marker), collapse = " OR "), ") AND ", min_length, ":", max_length, "[Sequence Length]", sep = "")
#    if(is.null(chunk_size)) {chunk_size=10000}
#  } else if (tolower(marker) %in% c("mitochondria", "mitochondrion")) {
#    message(paste0("Input marker is ", marker, ", Downloading full mitochondrial genomes"))
#    searchQ <- paste("(", x, "[ORGN])", " AND mitochondrion[filter] AND genome", sep = "")
#    if(is.null(chunk_size)) {chunk_size=1000}
#  } else if (tolower(marker) %in% c("genome", "mitochondrion")){
#    message(paste0("Input marker is ", marker, ", Downloading full genomes"))
#    searchQ <- paste("(", x, "[ORGN])", " AND complete genome[title]", sep = "")
#    if(is.null(chunk_size)) {chunk_size=1}
#  }
#  if(!quiet){message("Searching genbank with query:", searchQ)}
#
#  search_results <- rentrez::entrez_search(db = database, term = searchQ, retmax = 9999999, use_history = TRUE)
#  ids <- search_results$ids
#
#  #Convert search result GIDs to accession
#  if(length(ids) > 0 ){
#    accs <- gid_to_acc(ids, database = database, chunk_size = chunk_size,  multithread=multithread, progress = progress) %>%
#      stringr::str_remove(".[0-9]$")
#  } else{
#    stop("search returned no hits")
#  }
#  # Find any missing accessions
#  newsearch <- setdiff(accs, current)
#
#  # Download missing accessions
#
#  if (length(newsearch) > 0) {
#
#    if (!quiet) {message(paste0(length(newsearch), " sequences to be downloaded for: ", searchQ))}
#
#    chunks <- split(newsearch, ceiling(seq_along(newsearch)/chunk_size))
#
#    # setup multithreading
#    setup_multithread(multithread = multithread, quiet=quiet)
#
#    #Main function
#    furrr::future_map(chunks, function(x){
#      upload <- rentrez::entrez_post(db=database, id=x)
#      dl <- rentrez::entrez_fetch(db = database, web_history = upload, rettype = "gb", retmax = chunk_size)
#      gb <- biofiles::gbRecord(rcd = textConnection(dl))
#
#      # Hierarchial output
#      if (output == "standard") {
#        names <- paste0(biofiles::getAccession(gb), " ", biofiles::getDefinition(gb))
#      } else if (output == "h") {
#        lineage <- biofiles::getTaxonomy(gb) %>%
#          str_split_fixed(pattern = ";", n = Inf) %>%
#          trimws(which = "both") %>%
#          as_tibble() %>%
#          dplyr::mutate(Species = biofiles::getOrganism(gb)) %>%
#          dplyr::mutate(Genus = str_replace(V15, pattern = "[.]", replacement = "")) %>%
#          tidyr::unite("names", c("V1", "V2", "V4", "V6", "V10", "V14", "Genus", "Species"), sep = ";") %>%
#          dplyr::mutate(names = str_replace(names, pattern = " ", replacement = "_"))
#        names <- paste0(names(biofiles::getSequence(gb)), ";", lineage$names)
#
#        # Genbank taxID output
#      } else if (output == "gb") {
#        names <- cat_acctax(gb)
#      } else if (output == "binom") {
#        names <- paste0(biofiles::getAccession(gb), ";", biofiles::getOrganism(gb))
#      } else if (output == "gb-binom") {
#        names <- paste0(cat_acctax(gb), ";", biofiles::getOrganism(gb))
#      }
#      seqs <- biofiles::getSequence(gb) # Something failing here - getting much less than the proper amount
#      if(all(is.na(names(seqs)))) {
#        names(seqs) <- biofiles::getAccession(gb)
#      }
#      #Check if names match
#      names(seqs) <- names[names %>%
#                             stringr::str_remove(pattern=";.*$") %>%
#                             stringr::str_remove(pattern="\\|.*$") %>%
#                             stringr::str_remove(" .*$")
#                           %in% names(seqs)]
#      if (compress == TRUE) {
#        Biostrings::writeXStringSet(seqs, out_file, format = "fasta", compress = "gzip", width = 20000, append = TRUE)
#      } else if (compress == FALSE) {
#        Biostrings::writeXStringSet(seqs, out_file, format = "fasta", width = 20000, append = TRUE)
#      }
#      invisible(NULL)
#    }, .progress = progress)
#
#  }
#  #Close all workers
#  future::plan(future::sequential)
#
#  #Count number of downloaded sequences
#  if(file.exists(out_file)){
#    counter <- nrow(Biostrings::fasta.index(out_file))
#  } else {
#    counter <- 0
#  }
#
#  # Count total sequences that should have been downloaded
#  if(exists("search_results") && length(search_results$ids) > 1 ){
#    total_counter <- length(search_results$ids)
#  } else {
#    total_counter <- 0
#  }
#
#  res <- data.frame(
#    taxon = x,
#    seqs_total = total_counter,
#    seqs_downloaded = counter,
#    marker = marker,
#    database = database,
#    time = format(time, digits = 2)
#  )
#  return(res)
#}

#update_bold <- function(x, fasta, marker = "COI-5P", quiet = FALSE, output = "h", suffix="updates", compress = FALSE, force=FALSE,
#                        out_dir = NULL, db=NULL) {
#
#  # function setup
#  time <- Sys.time() # get time
#
#  if (!output %in% c("standard", "h", "binom", "gb", "bold", "gb-binom")) {
#    stop(paste0(output, " has to be one of: 'standard', 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
#  }
#
#  #Check if all inputs exists
#  if (!any(file.exists(fasta))) {
#    stop("Not all fasta files provided can be found, check the diferectory you have provided")
#  }
#
#  #Define directories
#  if (is.null(out_dir)) {
#    out_dir <- "bold"
#    if (!quiet) (message(paste0("No input out_dir given, saving output file to: ", out_dir)))
#  }
#  if (!dir.exists(out_dir)) {
#    dir.create(out_dir)
#  }
#
#  # Create output files
#  if (compress) {
#    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", marker, "_", suffix, ".fa.gz")
#  } else if (!compress) {
#    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", marker, "_", suffix, ".fa")
#  }
#
#  #Check if output exists
#  if (file.exists(out_file) && force==TRUE) {
#    file.remove(out_file)
#    cat("", file=out_file)
#  } else if (file.exists(out_file) && force==FALSE){
#    stop(paste0(out_file, " exists, set force = TRUE to overwrite"))
#  }
#
#  # Get NCBI taxonomy database if NCBI format outputs are desired
#  if (is.null(db) && output %in% c("gb", "gb-binom")) {
#    db <- get_ncbi_taxonomy(include_synonyms = TRUE, force=FALSE)
#  }
#
#  # Get accessions from existing fastas
#  current <- acc_from_fasta(fasta)
#
#  accs <-  bold::bold_specimens(taxon = x) %>%
#    dplyr::pull(sampleid)
#
#  newsearch <- setdiff(accs, current)
#
#  # Bold search
#  data <- bold::bold_seqspec(ids = newsearch, sepfasta = FALSE)
#
#  # Find differences
#
#  if (length(data) >0 && !methods::is(data, "logical")) {
#    data <- data %>%
#      dplyr::na_if("") %>%
#      dplyr::filter(stringr::str_detect(marker, markercode)) %>% # Remove all sequences for unwanted markers
#      dplyr::mutate(domain_name = "Eukaryota") %>%
#      dplyr::filter(!is.na(species_name)) %>%
#      dplyr::filter(!stringr::str_detect(nucleotides, pattern = "I")) #remove inosines
#
#    if (nrow(data) >0) {
#      if (output == "standard") {
#        data <- data %>%
#          dplyr::select(sampleid, species_name, markercode, genbank_accession, nucleotides) %>%
#          tidyr::unite("name", c("sampleid", "species_name", "markercode", "genbank_accession"), sep = "|")
#      } else if (output == "h") {
#        # Hierarchial output
#        data <- data %>%
#          dplyr::select(sampleid, domain_name, phylum_name, class_name,
#                        order_name, family_name, genus_name, species_name,nucleotides
#          ) %>%
#          tidyr::drop_na() %>%
#          tidyr::unite("name", c(
#            "sampleid", "domain_name",
#            "phylum_name", "class_name",
#            "order_name", "family_name",
#            "genus_name", "species_name"
#          ),
#          sep = ";"
#          ) %>%
#          dplyr::mutate(name = name %>%
#                          stringr::str_replace_all(pattern = " ", replacement = "_") %>%
#                          trimws(which = "both"))
#
#
#      } else if (output == "binom") {
#        # Binomial output
#        data <- data %>%
#          dplyr::select(sampleid, species_name, nucleotides) %>%
#          tidyr::drop_na() %>%
#          tidyr::unite("name", c("sampleid", "species_name"), sep = ";") %>%
#          dplyr::mutate(name = name %>%
#                          stringr::str_replace_all(pattern = " ", replacement = "_") %>%
#                          trimws(which = "both"))
#
#
#      } else if (output == "bold") {
#        # BOLD taxID output
#        data <- data %>%
#          dplyr::select(sampleid, species_taxID, nucleotides) %>%
#          tidyr::drop_na() %>%
#          tidyr::unite("name", c("sampleid", "species_taxID"), sep = ";") %>%
#          dplyr::mutate(name = name %>%
#                          stringr::str_replace_all(pattern = " ", replacement = "_") %>%
#                          trimws(which = "both"))
#
#      } else if (output == "gb") {
#        # Genbank taxID output
#        data <- data %>%
#          dplyr::select(sampleid, species_name, nucleotides) %>%
#          tidyr::drop_na() %>%
#          dplyr::mutate(tax_name = trimws(species_name, which = "both")) %>%
#          dplyr::left_join(db, by="tax_name") %>%
#          tidyr::unite("name", c("sampleid", "tax_id"), sep = "|")
#
#      } else if (output == "gb-binom") {
#        # genbank binomial output
#        data <- data %>%
#          dplyr::select(sampleid, species_name, nucleotides) %>%
#          tidyr::drop_na() %>%
#          dplyr::mutate(tax_name = trimws(species_name, which = "both")) %>%
#          dplyr::left_join(db, by="tax_name") %>%
#          tidyr::unite("name", c("sampleid", "tax_id"), sep = "|") %>%
#          tidyr::unite("name", c("name", "tax_name"), sep = ";")
#      }
#
#      seqs <- Biostrings::DNAStringSet(data$nucleotides)
#      names(seqs) <- data$name
#      if (compress == TRUE) {
#        Biostrings::writeXStringSet(seqs, out_file, format = "fasta", compress = "gzip", width = 20000)
#      } else if (compress == FALSE) {
#        Biostrings::writeXStringSet(seqs, out_file, format = "fasta", width = 20000)
#      }
#
#      # Done message
#      time <- Sys.time() - time
#      if (!quiet) (message(paste0("Downloaded ", length(seqs), " ", x, " Sequences from BOLD ", " in ", format(time, digits = 2))))
#    }
#  }
#  # Count number of downloaded sequences
#  if(file.exists(out_file)){
#    counter <- nrow(Biostrings::fasta.index(out_file))
#  } else {
#    counter <- 0
#  }
#
#  # Count total sequences that should have been downloaded
#  if(methods::is(data, "data.frame")){
#    total_counter <- nrow(data)
#  } else {
#    total_counter <- 0
#  }
#
#  res <- data.frame(
#    taxon = x,
#    seqs_total = total_counter,
#    seqs_downloaded = counter,
#    marker = marker,
#    database = "bold",
#    time = format(time, digits = 2)
#  )
#  return(res)
#}



# PHMM------------------------------------------------------------
#' Map to model
#'
#' @description This function alignes sequences to a Profile Hidden Markov Model (PHMM) using the Viterbi algorithm in order to retain only the target loci.
#' This function can also be used to extract smaller subregions out of longer sequences, for instance extracting the COI barcode region from mitochondrial genomes.
#' In order to reduce the number of sequences for alignment using the computationally expensive Viterbi algorithm, a rapid kmer distance screen is first conducted to remove any sequence too far diverged from the reference PHMM.
#' Similarly, all sequences that are more than twice the length of the reference PHMM model are broken into chunks of the same length as the model, and a rapid kmer screen conducted. The sequence is then subset to the most similar chunk and its two adjacent chunks prior to Viterbi alignment.
#'
#' @param x A DNAbin or DNAStringset object
#' @param model A Profile Hidden Markov model ("PHMM" object) generated by `aphid::derivePHMM` to align the sequences to.
#' An already derived model of COI can be loaded using `data("model", package="taxreturn")`
#' @param min_score The minimum specificity (log-odds score for the optimal alignment) between the query sequence and the PHMM model for the sequence to be retained.
#' see `?aphid::Viterbi` for more information about the alignment process.
#' @param min_length The minimum length of the match between the query and PHMM for a sequence to retained. Takes into account sequential matches, as well as any internal insertions or deletions below max_indel.
#' @param max_indel The maximum number of internal insertions or deletions within the sequence to allow.
#' @param max_gap The maximum number of gaps within the sequence to allow.
#' @param max_N The max number of ambiguous N bases allowed before a sequence is removed.
#' @param check_frame Whether sequences with insertions or deletions which arent in multiples of 3 should be removed from output. Useful for coding loci such as COI but should not be used for non-coding loci. Default is FALSE.
#' @param kmer_threshold the maximum kmer distance allowed from the reference model. If a sequence is further than this, it will be skipped from the slower Viterbi alignment. Default is 50% (0.5)
#' @param k integer giving the k-mer size used to generate the input matrix for k-means clustering. Default is k=5.
#' @param shave Whether bases that are outside (to the left or right) of the PHMM object should be removed from sequences in the output. Default is TRUE.
#' @param trim_ends Sometimes a trailing base can end up at the end of the alignment, separated by gaps. the trim_ends parameter checks up to n bases from each end of the alignment and if gaps are detected, any trailing bases will be removed.
#' @param extra How to handle insertions which were not part of the PHMM model. 'drop' will truncate all sequences to the shortest alignment length, while 'fill' will use gaps to pad all sequences out to the longest alignment length.
#' @param multithread Whether multithreading should be used, if TRUE the number of cores will be automatically detected (Maximum available cores - 1), or provide a numeric vector to manually set the number of cores to use. Default is FALSE (single thread)
#' @param quiet Whether progress should be printed to the console. Note that this will add additional runtime.
#' @param progress Whether a progress bar is displayed.
#'
#' @import stringr
#' @import furrr
#' @import future
#' @import dplyr
#' @importFrom ape as.DNAbin
#' @importFrom ape base.freq
#' @importFrom methods is
#'
#' @return
#' @export
#'
map_to_model <- function(x, model, min_score = 100, min_length = 1, max_indel = 9, max_gap = Inf, max_N = Inf,
                          check_frame = FALSE, kmer_threshold=0.5, k=5, shave = TRUE, trim_ends=FALSE, extra=NA,
                          multithread = FALSE, quiet = FALSE, progress = FALSE) {
  time <- Sys.time() # get time

  # Check inputs
  if (!methods::is(model, "PHMM")) {
    stop("Model needs to be a PHMM object")
  }
  if(!is.na(extra)){
    if(!extra %in% c("fill", "drop")) stop("extra must be 'fill', 'drop', or NA")
  }
  # Convert to DNAbin
  if (!methods::is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
    if (all(is.na(ape::base.freq(x)))) {stop("Error: Object is not coercible to DNAbin \n")}
  }
  if(!is.numeric(kmer_threshold) | !dplyr::between(kmer_threshold, 0,1)) { stop("Threshold must be a numeric between 0 and 1")}

  # Set up multithreading
  setup_multithread(multithread = multithread, quiet=quiet)

  # Ensure all sequences are above the kmer size
  no_gaps<- lapply(x, function(s) s[!(s %in% as.raw(c(2, 4, 240)))])
  x <- x[sapply(no_gaps, length) > k + 1]

  # Process long sequences & short sequences separately
  long_seqs <- x[lengths(x) > (model$size * 2)]
  short_seqs <- x[lengths(x) <= (model$size * 2)]

  # Conduct a quick kmer screen directly on the short sequences
  if(length(short_seqs) > 0){
    short_seqs <- kmer_screen(short_seqs, model, threshold=kmer_threshold, k = k, quiet= TRUE)
  }

  # Split the long sequences into chunks and conduct kmer screen on each chunk
  if(length(long_seqs) > 0){
    long_seqs <- long_seqs %>%
      furrr::future_map(subset_long_seq, model = model, split_length = model$size, threshold = kmer_threshold, k = k,
                        .progress = progress, .options = furrr::furrr_options(seed = TRUE))
    class(long_seqs) <- "DNAbin"
  }

  # Join short and long seqs together again
  x[names(long_seqs)] <- long_seqs
  x <- x[!names(x) %in% names(x)[!names(x) %in% c(names(short_seqs), names(long_seqs))]]

  #Apply filt_phmm to all sequences
  res <- furrr::future_map(x, filt_phmm, model=model, min_score=min_score, min_length=min_length, max_indel = max_indel,
                           shave=shave, trim_ends = trim_ends, check_frame=check_frame,
                           .progress = progress, .options = furrr::furrr_options(seed = TRUE))

  #Close all worker threads
  future::plan(future::sequential)

  # Filter sequences that didnt pass min_score or min_length
  discards <- sapply(res, is.null)

  if (sum(!discards) > 0) {
    if (!quiet) message("Retained ", sum(!discards), " sequences after alignment to PHMM \n")
    res <- res[!discards]
  } else {
    warning("None of the sequences met PHMM specificity criteria. Returning NULL \n")
    return(NULL)
  }

  # Remove sequences above max_N
  if(max_N < Inf){
    discards <- sapply(res, function(s) sum(s == as.raw(240))) > max_N
    res <- res[!discards]
    if (!quiet) message("Retained ", length(res), " sequences after applying ambiguity filter \n")
  }

  # Remove sequences above max_gap
  if(max_gap < Inf){
    discards <- sapply(res, function(s) sum(s == as.raw(c(4))))  > max_gap
    res <- res[!discards]
    if (!quiet) message("Retained ", length(res), " sequences after applying gap filter \n")
  }

  # Ensure outputs are consistent length
  size_range <- range(sapply(res, length))
  if(!(size_range[1]==size_range[2])){
    if(!is.na(extra) & extra=="fill"){
      res <- lapply(res, function(y) {
        z <- y[1:size_range[2]]
        z[z==00] <- as.raw(04)
        return(z)
      })
      message("Output sequences were of mixed lengths between ", size_range[1], "bp and ", size_range[2], "bp extra was filled")
    } else if(!is.na(extra) & extra=="drop"){
      res <- lapply(res, function(y) {y[1:size_range[1]]})
      message("Output sequences were of mixed lengths between ", size_range[1], "bp and ", size_range[2], "bp extra was dropped")
    } else if(is.na(extra)){
      warning("Output sequences are not all the same length due to insertion positions not present in the model, set extra='drop'truncate all sequences to the shortest alignment or extra='fill' to fill with gaps")
    }
  }
  class(res) <- "DNAbin"
  time <- Sys.time() - time
  if (!quiet) (message(paste0("finished in ", format(time, digits = 2))))
  return(res)
}

filt_phmm <- function(s, model, min_score = 100, min_length = 1, max_indel = Inf, shave = TRUE, check_frame=TRUE, trim_ends=FALSE, ...) {
  s <- s[!s %in% as.raw(c(2, 4))] #Remove gaps
  vit <- aphid::Viterbi(model, s, odds = TRUE, type = "semiglobal", cpp = TRUE, residues = "DNA", ...=...)

  # Stop early if below score or no matching states
  if (vit$score < min_score | !any(vit$path == 1)) {
    return(NULL)
  }

  # Calculate longest match length
  rle_path <- rle(vit$path)
  rle_matches <- rle_path$lengths
  names(rle_matches) <- rle_path$values

  # Only count alignments over 10 bases to prevent bad alignments
  rle_matches <- rle_matches[!names(rle_matches)=="1" | rle_matches > 10]

  # If there is more than one matching segment, see if the gaps between them are viable
  if(sum(names(rle_matches) == "1") > 1){
    all_matches <- rle_matches[names(rle_matches) == "1"]
    match_ranges <- range(which(names(rle_matches) == "1"))
    all_between <- seq(match_ranges[1], match_ranges[2], 1)
    gap_lengths <- rle_matches[all_between[!all_between %in% which(names(rle_matches) == "1")]]

    # Check if gap lengths is divisible by 3
    if(check_frame & !all((gap_lengths %% 3)==0, na.rm = TRUE)){
      return(NULL)
    }
    # If the numbers between the 1 elements are less than max_indel combine with previous. Starting from second element!
    for(m in 2:length(rle_matches)){
      # Check if the current element is not a 1, but there is a 1 ahead and behind
      if(all((!names(rle_matches[m]) == "1"), (names(rle_matches[m-1]) == "1"), (names(rle_matches[m+1]) == "1"), na.rm =TRUE)){
        #check if the current element is below max_indel
        if(rle_matches[m] <  max_indel){
          # if so, add the numbers
          names(rle_matches)[m] <- "1"
          rle_matches[m] <- (rle_matches[m] + rle_matches[m-1])
        }
      } else if(names(rle_matches[m]) == "1" & names(rle_matches[m-1]) == "1"){
        rle_matches[m] <- (rle_matches[m] + rle_matches[m-1])
      }
    }
  }
  # Get longest match
  longest_match <- max(rle_matches[names(rle_matches) == "1"])

  # Return null if below min_length
  if (longest_match < min_length) {
    return(NULL)
  }

  # if shave is true, shave to alignment 2 - Maybe here is where i need to put the rolling match call?
  if (shave & any(vit$path==2)) {
    last <- match(c(0, 1), rev(vit$path)) - 1
    last <- min(last[!is.na(last)])  # Getting the min, ignores deletions in the middle
    begin <- match(c(0, 1), vit$path)
    begin <- min(begin[!is.na(begin)])
    s <- s[begin: (length(s) - last)]
    vit$path <- vit$path[begin:(length(vit$path)-last)]
  }

  # Insert gaps to pad sequence out to the alignment length
  out <- as.raw(vit$path)
  out[out %in% c("01", "02")] <- s
  out[out==00] <- as.raw(04)

  # Trim any single bases on each end
  if (is.numeric(trim_ends)){
    if(out[trim_ends] == as.raw(04)){
      out[1:trim_ends] <- as.raw(04)
    }
    if(out[length(out)-(trim_ends-1)]== as.raw(04)){
      out[length(out)-(trim_ends-1):length(trim_ends)] <- as.raw(04)
    }
  }
  return(out)
}



# Codon filters ---------------------------------------------------------------
#' Get Reading frame of sequences
#'
#' @param x Sequences in DNAStringset or DNAbin format
#' @param genetic_code A genetic code for the Amino acid translation. set to 'SGC4' for Invertebrate mitochondrial or see all known codes at Biostrings::GENETIC_CODE_TABLE
#' @param tryrc Whether the reverse complemement should be evaluated if no frame without stop codons was found in the forward orientation.
#' @param resolve_draws How draws should be resolved when multiple possible frames produce sequences with no stop codons.
#' Options are "remove" to completely remove the sequence, or "majority" to pick the most common frame from the entire alignment.
#'
#' @return
#' @export
#'
#' @import stringr
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings translate
#' @importFrom Biostrings getGeneticCode
#' @importFrom methods is
#'
get_reading_frame <- function(x, genetic_code = NULL, tryrc=TRUE, resolve_draws="majority") {
  if(is.null(genetic_code)){
    stop("genetic_code must not be NULL, set to 'SGC4' for Invertebrate mitochondrial or see Biostrings::GENETIC_CODE_TABLE for other options")
  }
  # Convert to DNAStringSet
  if (methods::is(x, "DNAbin")) {
    x <- DNAbin2DNAstringset(x, remove_gaps=FALSE)
  }

  # Get reading frames for forward and reverse oriientations
  if(!tryrc){
    frames <- lapply(1:3, function(pos) XVector::subseq(x, start=pos))
  }else if (tryrc){
    frames <- c(lapply(1:3, function(pos) XVector::subseq(x, start=pos)),
                lapply(1:3, function(pos) XVector::subseq(Biostrings::reverseComplement(x), start=pos))
    )
  }
  #Translate all reading frames
  suppressWarnings(
    translated <- lapply(frames, Biostrings::translate, genetic.code = Biostrings::getGeneticCode(genetic_code), if.fuzzy.codon=c("solve", "X"))
  )

  #select the reading frames that contain 0 stop codons, or return NA
  reading_frame <- vector("integer", length=length(x))
  for (i in 1:length(x)){
    # Check forward frames first
    fvec <- c(stringr::str_count(as.character(translated[[1]][i]), "\\*"),
              stringr::str_count(as.character(translated[[2]][i]), "\\*"),
              stringr::str_count(as.character(translated[[3]][i]), "\\*"))
    if(sum(fvec==0)==1){
      # If only 1 appropriate frame for fwd direction
      reading_frame[i] <- which(fvec==0)
    } else if(sum(fvec==0)>1) {
      # If multiple appropriate frames for fwd direction
      reading_frame[i] <- 0
    } else if(sum(fvec==0)==0) {
      # If no appropriate frames for fwd direction
      # Try reverse direction or return NA
      if(tryrc){
        rvec <- c(stringr::str_count(as.character(translated[[4]][i]), "\\*"),
                  stringr::str_count(as.character(translated[[5]][i]), "\\*"),
                  stringr::str_count(as.character(translated[[6]][i]), "\\*"))
        if(sum(rvec==0)==1){
          #Check if only 1 appropriate frame for rev direction (return negative)
          reading_frame[i] <- -which(rvec==0)
        } else if(sum(rvec==0)>1) {
          #Check if multiple appropriate frames for rev direction
          reading_frame[i] <- 0
        }else if(sum(rvec==0)==0) {
          #Check if multiple appropriate frames for rev direction
          reading_frame[i] <- NA
        }
      } else {
        reading_frame[i] <- NA
      }
    }
  }
  # if resolve draws is majority, select the most common frame from across the whole dataset
  if (resolve_draws == "majority") {
    reading_frame[reading_frame==0] <- reading_frame[which.max(tabulate(reading_frame))]
  } else if (resolve_draws == "remove") {
    reading_frame[reading_frame==0] <- NA
  }
  return(reading_frame)
}


#' Filter sequences containing stop codons
#'
#' @param x Sequences in DNAStringset or DNAbin format
#' @param genetic_code A genetic code for the Amino acid translation. set to 'SGC4' for Invertebrate mitochondrial or see all known codes at Biostrings::GENETIC_CODE_TABLE
#' @param tryrc Whether the reverse complemement should be evaluated if no frame without stop codons was found in the forward orientation.
#' @param resolve_draws How draws should be resolved when multiple possible frames produce sequences with no stop codons.
#' Options are "remove" to completely remove the sequence, or "majority" to pick the most common frame from the entire alignment.
#'
#' @return
#' @export
#' @importFrom ape as.DNAbin
#' @importFrom Biostrings reverseComplement
#' @importFrom DECIPHER RemoveGaps
#' @importFrom methods is
#'
codon_filter <- function(x, genetic_code = NULL, tryrc=TRUE, resolve_draws="majority"){
  if(is.null(genetic_code)){
    stop("genetic_code must not be NULL, set to 'SGC4' for Invertebrate mitochondrial or see Biostrings::GENETIC_CODE_TABLE for other options")
  }
  # Convert to DNAStringSet
  if (methods::is(x, "DNAbin")) {
    x <- DNAbin2DNAstringset(x, remove_gaps=FALSE)
    format <- "DNAbin"
  } else if(methods::is(x, "DNAStringSet")){
    format <- "DNAStringSet"
  } else {
    stop("x must be a DNAbin or DNAStringSet")
  }

  #Get reading frames
  frames <- get_reading_frame(DECIPHER::RemoveGaps(x), genetic_code = genetic_code, tryrc = tryrc, resolve_draws = resolve_draws)

  # Check if any sequences need RC (negative reading frame)
  to_rc <- sign(frames)==-1
  to_rc[is.na(to_rc)] <- FALSE
  if (any(to_rc)){
    message(sum(to_rc), " Sequences reverse complemented as no forward match was found")
    x[to_rc] <- Biostrings::reverseComplement(x[to_rc])
  }

  if(format=="DNAbin"){
    out <- ape::as.DNAbin(x[!is.na(frames)])
  } else if(format=="DNAStringSet"){
    out <- x[!is.na(frames)]
  }
  message(paste0(length(x) - length(out), " Sequences containing stop codons removed"))
  return(out)
}

#' Codon entropy
#'
#' @param x Sequences in DNAStringset or DNAbin format
#' @param genetic_code A genetic code for the Amino acid translation. set to 'SGC4' for Invertebrate mitochondrial or see all known codes at Biostrings::GENETIC_CODE_TABLE
#' @param tryrc Whether the reverse complemement should be evaluated if no frame without stop codons was found in the forward orientation.
#' @param codon_filter Whether ` codon_filter` should be run first to remove sequences containing stop codons or frameshifts.
#' @param resolve_draws How draws should be resolved when multiple possible frames produce sequences with no stop codons.
#' Options are "remove" to completely remove the sequence, or "majority" to pick the most common frame from the entire alignment.
#' @param method the method employed to estimate entropy. see `?entropy::entropy` for more details
#'
#' @return
#' @export
#' @import purrr
#' @importFrom entropy entropy
#' @importFrom methods is
#'
#'
codon_entropy <- function(x, genetic_code = NULL, tryrc = TRUE, codon_filter = TRUE, resolve_draws = "majority", method = "ML") {
  if(is.null(genetic_code)){
    stop("genetic_code must not be NULL, set to 'SGC4' for Invertebrate mitochondrial or see Biostrings::GENETIC_CODE_TABLE")
  }
  if (methods::is(x, "DNAbin")) {
    x <- DNAbin2DNAstringset(x, remove_gaps=FALSE)
  }
  #Filter out sequences with stop codons
  if(codon_filter){
    x <- codon_filter(x, genetic_code = genetic_code, tryrc=tryrc)
  }

  #subset to the reading frame
  pos <- get_reading_frame(x, genetic_code = genetic_code, tryrc = tryrc, resolve_draws = resolve_draws)

  F_frames <-  as.character(XVector::subseq(x, start= pos))

  ent <- vector("list", length=length(F_frames))
  for (l in 1:length(F_frames)){
    ent[[l]] <- c(
      entropy::entropy(table(purrr::map_chr(seq(1, nchar(F_frames[l]), 3), function(i) substr(F_frames[l], i, i))), method=method),
      entropy::entropy(table(purrr::map_chr(seq(2, nchar(F_frames[l]), 3), function(i) substr(F_frames[l], i, i))), method=method),
      entropy::entropy(table(purrr::map_chr(seq(3, nchar(F_frames[l]), 3), function(i) substr(F_frames[l], i, i))), method=method)
    )
    names(ent[[l]]) <- c("pos1", "pos2", "pos3")
  }
  names(ent) <- names(x)
  return(ent)
}


# kmer screening -------------------------------------------------------------
#' Kmer screening function
#'
#' @param seqs A DNAbin object
#' @param model A Profile Hidden Markov model ("PHMM" object) generated by `aphid::derivePHMM` to screen the sequences against.
#' An already derived model of COI can be loaded using `data("model", package="taxreturn")`
#' @param threshold The maximum kmer distance allowed from the reference model before a sequence is discarded.
#' @param k The k-mer size to be used for calculating the distance matrix, parsed to kmer::mbed. Note that high values of k may be slow to compute and use a lot of memory due to the large numbers of calculations required, particularly when the residue alphabet is also large.
#' @param quiet Whether progress should be printed to the console.
#'
#' @return
#' @export
#' @importFrom kmer mbed
#' @importFrom aphid generate
#'
kmer_screen <- function(seqs, model, threshold = 0.3, k = 5, quiet=FALSE){
  # Check inputs
  if(!class(model) == "PHMM") {
    stop(paste0("Model must be of class 'PHMM', but it is of class '",class(model),"'."))
  }
  if(!class(seqs) == "DNAbin") {
    stop(paste0("Seqs must be of class 'DNAbin', but it is of class '",class(seqs),"'."))
  }

  # Generate a reference sequence from the model to align against
  seed_seq <- char2DNAbin(paste(aphid::generate(model, model$size, random=FALSE), collapse=""))
  names(seed_seq) <- "SEED"

  # Align against seed sequence using mbed
  mbed_dist <- kmer::mbed(c(seed_seq, seqs), seeds = 1, k = k, residues = NULL, gap = "-",
                          counts = FALSE)[,]

  # Return only sequences below kmer distance threshold
  seqs_to_keep <- names(mbed_dist)[mbed_dist < threshold]
  seqs_to_keep <- seqs_to_keep[seqs_to_keep %in% names(seqs)]

  out <- seqs[seqs_to_keep]
  if(length(out) > 0){
    if(!quiet){ message(length(seqs[!names(seqs) %in% names(out)]), " sequences >", threshold, " ", k,"-mer distance removed")}
    return(out)
  } else{
    if(!quiet){ warning("All sequences are >", threshold, " ", k,"-mer distance and were removed")}
    return(NULL)
  }
  return()
}

#' Closest kmer distance
#' This function checks the kmer distance between all the sequences in a DNAbin and a reference PHMM model
#'
#' @param seqs A DNAbin object
#' @param model A Profile Hidden Markov model ("PHMM" object) generated by `aphid::derivePHMM` to screen the sequences against.
#' An already derived model of COI can be loaded using `data("model", package="taxreturn")`
#' @param threshold The maximum kmer distance allowed from the reference model before a sequence is discarded.
#' @param k The k-mer size to be used for calculating the distance matrix, parsed to kmer::mbed. Note that high values of k may be slow to compute and use a lot of memory due to the large numbers of calculations required, particularly when the residue alphabet is also large.
#' @param quiet Whether progress should be printed to the console.
#'
#' @return
#' @export
#'
#' @importFrom kmer mbed
#' @importFrom aphid generate
#'
closest_seq <- function(seqs, model, threshold, k=5, quiet=FALSE){
  # Check inputs
  if(!class(model) == "PHMM") {
    stop("Model must be of class PHMM")
  }
  if(!class(seqs) == "DNAbin") {
    stop("Seqs must be of class DNAbin")
  }

  # Generate a reference sequence from the model to align against
  seed_seq <- char2DNAbin(paste(aphid::generate(model, model$size), collapse=""))
  names(seed_seq) <- "SEED"

  # Align against seed sequence using mbed
  mbed_dist <- kmer::mbed(c(seed_seq, seqs), seeds = 1, k = k, residues = NULL,
                          counts = FALSE)[,]
  if(any(mbed_dist < threshold)){
    out <- unname(which.min(mbed_dist[!names(mbed_dist) == "SEED"]))
  } else (out <- NULL)
  return(out)
}

chunk_seqs <- function(x, split_length){
  n <- round(length(x)/split_length)
  chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
  if(n >1 ){
    splits <- chunk(x, n)
  } else (splits <- x)

  #Check no chunks are just all gaps
  no_gaps<- lapply(splits, function(s) s[!(s %in% as.raw(c(2, 4, 240)))])
  splits <- splits[sapply(no_gaps, length) > 0]

  class(splits) <- "DNAbin"
  return(splits)
}

subset_long_seq <- function(x, model, split_length, threshold=0.3, k=5, quiet=FALSE){
  # Split sequences into chunks
  splits <- chunk_seqs(x, split_length=split_length)
  if(!is.list(splits))(return(x))

  # Check all chunks are above value of k
  no_gaps<- lapply(splits, function(s) s[!(s %in% as.raw(c(2, 4, 240)))])
  splits <- splits[sapply(no_gaps, length) > (k+1)]

  # find chunk with closest kmer distance
  min_split <- closest_seq(splits, model = model, threshold = threshold, k = k, quiet = quiet)

  # Return that chunk and its adjacent neighbours - within bounds of splits
  chunks_to_keep <- c(min_split-1,min_split, min_split+1)
  chunks_to_keep <- chunks_to_keep[(chunks_to_keep > 0 & chunks_to_keep <= length(splits))]

  # Return concatenated subsequence
  out <- unlist(splits[chunks_to_keep])
  return(out)
}

# Prune group sizes -------------------------------------------------------
#' Prune group sizes
#'
#' @param x A DNAbin or DNAStringset object
#' @param max_group_size The maximum number of sequences with the same taxonomic annotation to keep
#' @param dedup Whether sequences with identical taxonomic name and nucleotide bases sequences should be discarded first
#' @param discardby How sequences from groups with size above max_group_size should be discarded.
#' Options include "length" (Default) which will discard sequences from smallest to largest until the group is below max_group_size,
#' "random" which will randomly pick sequences to discard until the group is below max_group_size.
#' @param prefer A vector of sequence names that will be preferred when subsampling groups when discardby=random,
#' or prefered when breaking ties in sequences of the same length when discardby=length. For instance high quality in-house sequences.
#' @param quiet Whether progress should be printed to the console.
#' @return
#' @export
#'
#' @import dplyr
#' @import stringr
#' @importFrom tibble as_tibble
#' @importFrom tidyr separate
#' @importFrom openssl md5
#' @importFrom ape as.DNAbin
#' @importFrom ape base.freq
#' @importFrom methods is
#'
prune_groups <- function(x, max_group_size = 5, dedup = TRUE, discardby = "length", prefer=NULL, quiet = FALSE) {
  # Convert to DNAbin
  if (!methods::is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
    if (all(is.na(ape::base.freq(x)))) {stop("Error: Object is not coercible to DNAbin \n")}
  }
  if (!is.null(prefer) & !is.character(prefer)){
    stop("prefer must be either NULL or a vector of sequence names to prefer")
  }

  # Remove duplicate sequences
  if (dedup) {
    dup <- length(x)
    taxids <- names(x) %>%
      stringr::str_remove("^.*;")

    ## Consider taxonomic name and sequence identity in deduplication
    hashes <- purrr::map2(x, taxids, ~{
      openssl::md5(paste(c(.y, as.vector(.x)), collapse=""))
    }) %>%
      unlist()

    # get list of all duplicated hashes
    dupes <- unique(hashes[duplicated(hashes)])
    remove <- logical(length(x))
    for (i in 1:length(dupes)) {
      index <- which(hashes %in% dupes[i])
      # Check if in prefered - If so keep first element in the prefered
      if (is.character(prefer) & any(names(x)[index] %in% prefer)){
        keep <- index[names(x)[index] %in% prefer][1]
      } else {
        # otherwise just pick first element
        keep <- index[1]
      }
      remove[index[!index %in% keep]] <- TRUE
    }
    x <- x[!remove]
    if (!quiet) cat(paste0((dup - length(x)), " duplicate sequences removed \n"))
  }

  # Remove sequences from groups where more species names than max_group_size
  groups <- names(x) %>%
    stringr::str_split_fixed(";", n = 2) %>%
    tibble::as_tibble(.name_repair = ~ c("acc", "taxon")) %>%
    dplyr::pull(taxon) %>%
    openssl::md5()
  groupCounts <- table(groups) # Count number of seqs per group
  u_groups <- names(groupCounts) # Get unique groups

  remove <- logical(length(x))
  # Random discarding
  if (discardby == "random") {
    for (i in which(groupCounts > max_group_size)) {
      index <- which(groups == u_groups[i])
      # Deal with prefered
      if (is.character(prefer) & any(names(x)[index] %in% prefer)){
        in_pref <- index[names(x)[index] %in% prefer]
        # Check if can sample just from prefered, or need a mix of the two
        if(in_pref >= max_group_size){
          keep <- sample(
            x = in_pref,
            size = max_group_size,
            replace = FALSE
          )
        } else{
          to_sample <- max_group_size - length(in_pref)
          keep <- c(in_pref,
                    sample(
                      x = index[-in_pref],
                      size = to_sample,
                      replace = FALSE
                    ))
        }
      } else {
        # otherwise just take random sample
        keep <- sample(
          x = index,
          size = max_group_size,
          replace = FALSE
        )
      }
      remove[index[-keep]] <- TRUE
    }
    # length discarding
  } else if (discardby == "length") {
    for (i in which(groupCounts > max_group_size)) {
      index <- which(groups == u_groups[i])
      rem <- sapply(x[index], function(s) length(s) - sum(s == as.raw(c(4)))) # Get lengths
      names(rem) <- index
      rem <- sort(rem, decreasing = TRUE)

      # deal with prefered
      if (is.character(prefer) & any(names(x)[index] %in% prefer)){
        in_pref <- index[names(x)[index] %in% prefer]
        non_pref <- index[!names(x)[index] %in% prefer]

        # Get the minimum length of sequences that are being dropped
        min_len <- min(rem[1:max_group_size])

        # Check if there are mixed prefered and non prefered at that length
        mixed_pref <- rem[names(rem) %in% in_pref] == min_len
        mixed_pref <- names(mixed_pref)[mixed_pref]
        mixed_non_pref <- rem[names(rem) %in% non_pref] == min_len
        mixed_non_pref <- names(mixed_non_pref)[mixed_non_pref]

        if((length(mixed_pref) > 0) & (length(mixed_pref) > 0)){
          # re-sort vector to prefer
          sample_from <- as.integer(c(names(rem)[rem > min_len], mixed_pref, mixed_non_pref, names(rem)[rem < min_len]))
          keep <- as.integer(sample_from[1:max_group_size])
        } else {
          # otherwise just take the longest
          keep <- as.integer(names(rem[1:max_group_size]))
        }
      } else{
        keep <- as.integer(names(rem[1:max_group_size]))
      }
      remove[index[!index %in% keep]] <- TRUE
    }
  }
  x <- x[!remove]
  if (!quiet) cat(paste0(sum(remove), " sequences pruned from over-represented groups"))
  return(x)
}

# Mixed clusters -----------------------------------------------------
#' Get mixed clusters
#'
#' @description Cluster sequences at a certain taxonomic similarity, and find clusters that contain mixed taxonomic names,
#' @description Note, it is recommended to set a unique seed using set.seed()
#' @param x	 A DNAbin list object whose names include NCBItaxonomic identification numbers.
#' @param db A taxonomic database from `get_ncbi_taxonomy` or `get_ott_lineage`
#' @param rank The taxonomic rank to check clusters at, accepts a character such as "order", or vector of characters such as c("species", "genus").
#' If "all", the clusters will be checked at all taxonomic ranks available.
#' @param threshold numeric between 0 and 1 giving the OTU identity cutoff for clustering. Defaults to 0.97.
#' @param rngseed (Optional) A single integer value passed to set.seed, which is used to fix a seed for reproducibly random number generation for the kmeans clustering.
#' If set to FALSE, then no fiddling with the RNG seed is performed, and it is up to the user to appropriately call set.seed beforehand to achieve reproducible results.
#' @param confidence The minimum confidence value for a mixed cluster to be flagged. For example, if confidence = 0.8 (the default value)
#'  a cluster will only be flagged if the taxonomy of a sequence within the cluster differs from at least four other independent sequences in its cluster.
#'  @param nstart how many random sets should be chosen for `kmeans`, It is recommended to set the value of nstart to at least 20.
#'  While this can increase computation time, it can improve clustering accuracy considerably.
#' @param return What type of data about the data should be returned. Options include:
#' Consensus - The consensus taxonomy for each cluster and associated confidence level
#' All - Return all taxa in mixed clusters and their sequence accession numbers
#' Count - Return counts of all taxa within each cluster
#' @param k integer giving the k-mer size used to generate the input matrix for k-means clustering.
#' @param quiet logical indicating whether progress should be printed to the console.
#' @param ... further arguments to pass to kmer::otu.
#'
#' @return
#' @export
#' @import dplyr
#' @importFrom kmer otu
#' @importFrom tibble rownames_to_column
#' @importFrom methods as
#' @examples
#' \dontrun{
#' seqs <- ape::read.FASTA("test.fa.gz")
#'
#' # NCBI taxonomy
#' mixed <- get_mixed_clusters(seqs, db, rank="species", threshold=0.99, confidence=0.8, quiet=FALSE)
#'
#' # OTT taxonomy
#' seqs <- map_to_ott(
#' seqs, dir="ott3.2", from="ncbi",
#' resolve_synonyms=TRUE, filter_bads=TRUE, remove_na = TRUE, quiet=FALSE
#' )
#'
#' mixed <- get_mixed_clusters(
#' seqs, db, rank="species",
#' threshold=0.99, confidence=0.6, quiet=FALSE
#' )
#' }
#'
get_mixed_clusters <- function (x, db, rank = "order", threshold = 0.97, rngseed = FALSE, confidence = 0.8, return = "consensus", k=5, quiet = FALSE, ...) {
  if(missing(x)) {stop("Error: x is required")}

  #Check inputs
  if(!is.numeric(threshold) | !dplyr::between(threshold, 0,1)) { stop("Threshold must be a numeric between 0 and 1")}
  if(!is.numeric(confidence) | !dplyr::between(confidence, 0,1)) { stop("Confidence must be a numeric between 0 and 1")}
  rank <- tolower(rank)
  return <- tolower(return)
  if(!return %in% c("consensus", "all", "counts")){stop("Return must be one of: 'consensus', 'all', or 'counts'")}
  if (methods::as(rngseed, "logical")) {
    set.seed(rngseed)
    if (!quiet) {
      message("`set.seed(", rngseed, ")` was used to initialize repeatable random subsampling.")
      message("Please record this for your records so others can reproduce.")
      message("Try `set.seed(", rngseed, "); .Random.seed` for the full vector",
              sep = "")
    }
  } else if (!quiet) {message("You set `rngseed` to FALSE. Make sure you've set & recorded\n",
                              " the random seed of your session for reproducibility.\n",
                              "See `?set.seed`\n")}

  #Get lineage
  if(attr(db, "type")  == "ncbi"){
    source <- "ncbi"
    lineage <-  get_ncbi_lineage(x = x, db = db)
  } else if(attr(db, "type")  == "OTT"){
    source <- "OTT"
    lineage <-  get_ott_lineage(x = x, db = db)
  } else (stop("db type is not supported"))

  # Cluster OTUS
  if (is.null(attr(x, "OTU"))) {
    if (!quiet) {cat(paste0("Clustering OTUs at ", (threshold*100), "%  similarity \n"))}
    otus <- kmer::otu(x, k=k , threshold = threshold, ...)
  } else {
    if (!quiet) {cat("Obtaining OTU membership from input object\n")}
    otus <- attr(x, "OTU")
    stopifnot(length(x) == length(otus))
  }
  if(length(unique(otus))==1) {
    warning("Only one unique cluster")
    return(NULL)
  }
  if (!quiet) {cat("Comparing lineage metadata within OTUs\n")}

  # Get mixed clusters
  find_mixed <- function(y, return) {
    # Ensure no duplicated accessions+names to bias confidence
    hashes <- paste0(gsub("\\|.*$", "\\1", names(y)), y)
    yu <- y[!duplicated(hashes)]
    if (length(unique(yu)) < 2) {
      return(NULL)
    }
    # Tabulate taxon names
    tab <- sort(table(yu), decreasing = TRUE)

    if(return == "consensus"){
      consensus <- names(tab)[1]
      consensus_taxid <- gsub("^.*\\|", "\\1", names(y)[y==consensus][1])
      mixed <- y != consensus #potential misannotated
      mixedu <- yu != consensus
      nu <- length(mixedu)

      #Check if there is a clear consensus
      if (tab[1] == tab[2]){
        consensus <- NA
        consensus_taxid <- NA
      }
      res <- data.frame(listed = y[mixed],
                        consensus = rep(consensus, sum(mixed)),
                        consensus_taxid = rep(consensus_taxid,sum(mixed)),
                        confidence = sum(!mixedu)/nu, cluster_size = nu,
                        stringsAsFactors = FALSE)  %>%
        tibble::rownames_to_column("Acc")
    } else if (return=="counts"){
      res <- data.frame(tax_name = names(tab),
                        count =as.numeric(tab), stringsAsFactors = FALSE)
    } else if (return=="all"){
      res <- data.frame(Acc =  gsub("\\|.*$", "\\1", names(y)),
                        tax_id = gsub("^.*\\|", "\\1", names(y)),
                        listed = as.character(y), stringsAsFactors = FALSE)
    }

    return(res)
  }

  # Loop over rank
  results <- vector("list", length=length(rank))
  for (i in 1:length(rank)){
    lins <- lineage %>%
      dplyr::select(!!rank[i]) %>%
      dplyr::pull(!!rank[i])
    if(length(lins[is.na(lins)]) >0 & !quiet){warning("ignoring ", length(lins[is.na(lins)]), " sequence/s with no taxonomic data for ", rank[i])}
    names(lins) <- lineage$Acc
    f <- as.factor(otus[!is.na(lins)])
    lins <- lins[!is.na(lins)]

    splitlist <- split(lins, f)
    splitlist <- splitlist[tabulate(f) > 2]

    mixedtab <- lapply(splitlist, find_mixed, return)
    mixedtab <- mixedtab[!vapply(mixedtab, is.null, logical(1))]
    if (length(mixedtab) == 0) {
      if (!quiet) {cat("No mixed clusters at", rank[i],   "rank \n")}
      results[[i]] <-  as.data.frame(NULL)
    } else if (length(mixedtab) > 0){

      mixedtab <- dplyr::bind_rows(mixedtab, .id="cluster")
      if(rank[i] == "species"){
        if(return == "consensus"){
          #Return binomials if species rank is selected
          mixedtab <- mixedtab %>%
            dplyr::left_join(lineage %>% dplyr::select(Acc, genus), by="Acc") %>%
            dplyr::mutate(consensus_taxid = as.numeric(na_if(consensus_taxid, "NA"))) %>%
            dplyr::left_join(db %>% dplyr::select(consensus_taxid = tax_id, consensus_genus = genus), by="consensus_taxid") %>%
            dplyr::mutate(listed = paste0(genus, " ", listed),
                          consensus  = paste0(consensus_genus," ", consensus)) %>%
            dplyr::select(-genus, -consensus_genus)
        } else if(return == "all"){
          mixedtab <- mixedtab %>%
            dplyr::left_join(lineage %>%
                               dplyr::select(Acc, genus) %>%
                               tidyr::separate(Acc, into=c("Acc", "tax_id"), sep="\\|"), by=c("Acc", "tax_id")) %>%
            dplyr::mutate(listed = paste0(genus, " ", listed)) %>%
            dplyr::select(-genus)
        }
      }
      if(return == "consensus"){
        mixedtab <- mixedtab[mixedtab$confidence >= confidence, ]
      }
      if (nrow(mixedtab) == 0) {
        if (!quiet) {cat("No mixed clusters at", rank[i],   "rank \n")}
        results[[i]] <-  as.data.frame(NULL)
      } else if(nrow(mixedtab) > 0 ) {
        if(return == "consensus"){
          mixedtab <- mixedtab[order(mixedtab$confidence, decreasing = TRUE), ]
        }
        if (!quiet) {cat("identified", length(unique(mixedtab$cluster)), "mixed clusters at", rank[i],   "rank \n")}
        results[[i]] <-  mixedtab %>%
          dplyr::mutate(rank = rank[i],
                        threshold = threshold) %>%
          dplyr::mutate_if(is.factor, as.character)

      }
    }
  }

  out <- dplyr::bind_rows(results)
  if (nrow(out)==0) {
    return(NULL)
  } else (return(out))
}




# Alignment utilities -----------------------------------------------------

#' Get primer binding position
#'
#' @param primer A character string, DNAbin object or an object coercible to DNAbin
#' @param model A profile hidden Markov model (a "PHMM" object) generated by the aphid R package to align the sequences to.
#' @param tryrc Whether the reverse complement should also be aligned. The highest scoring complement is chosen.
#' @param quiet Whether progress should be printed to the console.
#' @param min_score Minimum score for the viterbi alignment.
#' @param ... aditional arguments to be passed to "Viterbi"
#'
#' @return
#' @export
#' @import stringr
#' @importFrom ape complement
#' @importFrom aphid Viterbi
#'
get_binding_position <- function (primer, model, tryrc = TRUE, quiet = FALSE, min_score = 10, ...) {

  input <- primer
  if (!inherits(model, "PHMM")) { stop("Error: model must be a PHMM object")}

  if (!is.null(primer)) {
    if (!inherits(primer, "DNAbin")) {
      if (mode(primer) == "character") {
        if (nchar(primer[1]) == 1) {primer <- paste0(primer, collapse = "")}
        if(stringr::str_detect(primer, "I")) {message(paste0("Warning: Inosine (I) bases detected in primer ", input," these will be converted to N!"))}
        primer <- char2DNAbin(primer)
      }
      else {
        if (!inherits(primer, "PHMM"))
          stop("Invalid primer(s)\n")
      }
    }
  }

  up <- primer[!primer %in% as.raw(c(2, 4))]
  vitF <- aphid::Viterbi(model, up, odds = TRUE, type = "semiglobal", cpp = TRUE, residues = "DNA", ...=...)

  # Calculate longest match length
  rle_path <- rle(vitF$path)
  rle_matches <- rle_path$lengths
  names(rle_matches) <- rle_path$values
  all_matches <- rle_matches[names(rle_matches) == "1"]

  # Check for primers split into multiple matches
  if(length(all_matches) > 1){
    longest_match <- max(all_matches)

    if(which.max(all_matches) == 1){
      up[[1]] <- up[[1]][1:max(all_matches)]
    } else {
      up[[1]] <- up[[1]][all_matches[which.max(all_matches)-1]:all_matches[which.max(all_matches)]]
    }
    # Redo the viterbi alignment with just the subset primer
    vitF <- aphid::Viterbi(model, up, odds = TRUE, type = "semiglobal", cpp = TRUE, residues = "DNA", ...=...)
  }

  if (tryrc == TRUE) {
    down <- ape::complement(primer)
    down <- down[!down %in% as.raw(c(2, 4))]
    vitR <- aphid::Viterbi(model, down, odds = TRUE, type = "semiglobal", cpp = TRUE, residues = "DNA", ...=...)
    match_type <- NA_character_

    # Calculate longest match length
    rle_path <- rle(vitR$path)
    rle_matches <- rle_path$lengths
    names(rle_matches) <- rle_path$values
    all_matches <- rle_matches[names(rle_matches) == "1"]

    # Check for primers split into multiple matches
    if(length(all_matches) > 1){
      longest_match <- max(all_matches)

      if(which.max(all_matches) == 1){
        down[[1]] <- down[[1]][1:max(all_matches)]
      } else {
        down[[1]] <- down[[1]][all_matches[which.max(all_matches)-1]:all_matches[which.max(all_matches)]]
      }
      # Redo the viterbi alignment with just the subset primer
      vitR <- aphid::Viterbi(model, down, odds = TRUE, type = "semiglobal", cpp = TRUE, residues = "DNA", ...=...)
    }
    if (vitF$score > vitR$score & vitF$score > min_score) {
      match_type <- "forward"
      path <- vitF$path
      score <- vitF$score
    } else if (vitF$score < vitR$score & vitR$score > min_score) {
      match_type <- "reverse"
      path <- vitR$path
      score <- vitR$score
    } else if (vitF$score < min_score & vitR$score < min_score) {
      score <- max(vitF$score, vitR$score)
      out <- data.frame(primer = input, start = NA, end = NA, score=score)
      warning("Both complements of primer were below min_score")
      return(out)
    }
  } else if(tryrc == FALSE & vitF$score > min_score) {
    path <- vitF$path
    score <- vitF$score
    match_type <- "forward"
  } else {
    out <- data.frame(primer = input, start = NA, end = NA, score=vitF$score)
    warning("Forward complement of primer was below min_score")
    return(out)
  }

  matchF <- match(1, path)
  matchR <- (length(path) - (match(1, rev(path)) - 1))
  if (!quiet) {
    if(match_type == "forward"){
      message(paste0("Forward complement matched alignment at positions ", matchF, ":",matchR))
    } else if(match_type == "reverse"){
      message(paste0("Reverse complement matched alignment at positions ", matchF, ":",matchR))
    } else {
      message("Both complements did not match alignment")
    }
  }

  if ((matchR - (matchF - 1)) == length(primer[[1]])) {
    out <- data.frame(primer = input, start = matchF, end = matchR, score=score)
  }  else if ((matchR - (matchF - 1)) > length(primer[[1]])) {
    warning("Binding positions are larger than the primer length")
    out <- data.frame(primer = input, start = matchF, end = matchR, score=score)
  }  else if ((matchR - (matchF - 1)) < length(primer[[1]])) {
    warning("Binding positions are less than the primer length")
    out <- data.frame(primer = input, start = matchF, end = matchR, score=score)
  }
  return(out)
}

#' Get subalignment
#'
#' @description Aligns a DNABin to a reference PHMM model, and returns the optimal path
#' @param x A DNAbin object or an object coercible to DNAbin
#' @param model A profile hidden Markov model (a "PHMM" object) generated by the aphid R package to align the sequences to.
#' @param tryrc Whether the reverse complement should also be aligned
#' @param quiet Whether progress should be printed to the console.
#' @param check_indels Check that indels are multiples of 3, recommended for coding sequences such as COI
#' @param min_score Minimum score for the viterbi alignment
#' @param omit_endgaps Should gap characters at each end of the sequences be ignored when deriving the transition probabilities of the model? Defaults to FALSE. Set to TRUE if x is not a strict global alignment (i.e. if the alignment contains partial sequences with missing sections represented with gap characters).
#' @param multithread Whether multithreading should be used, if TRUE the number of cores will be automatically detected, or provided a numeric vector to manually set the number of cores to use
#' @param ... aditional arguments to be passed to "Viterbi"
#'
#' @return
#' @export
#' @import future
#' @importFrom ape as.DNAbin
#' @importFrom ape base.freq
#' @importFrom ape complement
#' @importFrom aphid derivePHMM
#' @importFrom aphid Viterbi
#'
#'
get_subalignment <- function(x, model, tryrc=FALSE, quiet=FALSE, check_indels=TRUE, min_score=10, omit_endgaps	= FALSE, multithread=FALSE, ...) {

  # Ensure x is a DNAbin
  if (!inherits(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
    if (all(is.na(ape::base.freq(x)))) {
      stop("Error: Object is not coercible to DNAbin \n")
    }
  }

  # setup multithreading
  setup_multithread(multithread = multithread, quiet=quiet)

  up <- aphid::derivePHMM(x, cores = cores, omit.endgaps = omit_endgaps, ...=...)
  vitF <- aphid::Viterbi(model, up, odds = TRUE, type = "semiglobal", cpp = TRUE, residues = "DNA")

  # Derive PHMM
  if(tryrc == TRUE) {
    if(!quiet){message("Deriving PHMM for reverse complement")}
    down <- aphid::derivePHMM(ape::complement(x), cores = cores, omit.endgaps = omit_endgaps, ...=...)
    vitR <- aphid::Viterbi(model, down, odds = TRUE, type = "semiglobal", cpp = TRUE, residues = "DNA")

    if (vitF$score > vitR$score && vitF$score > min_score) {
      if(!quiet){message("Forward complement matched alignment")}
      path <- vitF$path

    } else if (vitF$score < vitR$score && vitR$score > min_score) {
      if(!quiet){message("Reverse complement matched alignment")}
      path <- vitR$path

    } else if( vitF$score && vitR$score < min_score ){
      return(NULL)
      stop("Both complements of primer were below min_score")

    }
  } else if(tryrc == FALSE && vitF$score > min_score) {
    path <- vitF$path
  } else {Stop("Forward complement of primer was below min_score")}
  return(path)
}


#' Pad alignment
#'
#' @description Aligns a DNABin to a reference PHMM model, and pads any gaps between the query and reference
#' @param x A DNAbin object or an object coercible to DNAbin
#' @param model A profile hidden Markov model (a "PHMM" object) generated by the aphid R package to align the sequences to.
#' @param pad The character used to pad the gaps
#' @param tryrc Whether the reverse complement should also be aligned
#' @param quiet Whether progress should be printed to the console.
#' @param check_indels Check that indels are multiples of 3, recommended for coding sequences such as COI
#' @param min_score Minimum score for the viterbi alignment
#' @param omit_endgaps Should gap characters at each end of the sequences be ignored when deriving the transition probabilities of the model? Defaults to FALSE. Set to TRUE if x is not a strict global alignment (i.e. if the alignment contains partial sequences with missing sections represented with gap characters).
#' @param multithread Whether multithreading should be used, if TRUE the number of cores will be automatically detected, or provided a numeric vector to manually set the number of cores to use
#' @param ... aditional arguments to be passed to "Viterbi"
#'
#' @return
#' @export
#' @import stringr
#'
#'
pad_alignment <- function(x, model, pad = "-", tryrc = FALSE, check_indels = TRUE, min_score = 10, omit_endgaps	= FALSE, multithread = FALSE,  quiet = FALSE, ...){

  path <- get_subalignment(x = x, model = model, tryrc = tryrc, check_indels = check_indels,
                           min_score = min_score, omit_endgaps=omit_endgaps, multithread = multithread, quiet = quiet, ...=...)

  # Find start, stop, and indels
  matchF <- match(2, path)
  matchR <- (length(path) - (match(2, rev(path))-1))
  indels <- which(path == 0, arr.ind=FALSE)
  # potentially could have indels as 1's as well, due to potential for indels to be recorded in PhMM?
  # Do another check for which(path ==1, arr.ind=FALSE), and make sure they are more than matchF and less than matchR to be recorded as indels

  # Pad Left
  if (matchF > 1){
    left_pad <- which(path == 1, arr.ind=FALSE)
    left_pad <- left_pad[which(left_pad < matchF)]
  } else(left_pad <- NULL)

  # Detect indels
  if (length(indels) > 1) {
    # Detect multiple indels
    splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
    split <- splitAt(indels, which(diff(indels) > 1, arr.ind=FALSE)+1)

    # Confirm all indels are 1 codon deletions
    if (check_indels == TRUE && !all(sapply(split, length) %in% seq(3,12,3))) {
      stop("ERROR: Indels are not in multiples of 3!")
    }
  } else (indels <- NULL)

  # Pad right
  if (matchR < length(path)){
    right_pad <- which(path == 1, arr.ind = FALSE)
    right_pad <- right_pad[which(right_pad > matchR)]
  } else (right_pad <- NULL)

  insert <- c(left_pad, indels+length(left_pad), right_pad+length(left_pad)+length(indels))

  x <- as.character(x)
  insert_at <- function(x, index) {
    x <-  paste0(x, collapse = "")
    for(i in 1:length(index)) {
      stringr::str_sub(x, start = index[i], end = index[i]-1) <- pad
    }
    x <- stringr::str_to_upper(x)
    return(x)
  }
  out <- sapply(x, insert_at, insert)
  out <- char2DNAbin(out)
  return(out)
}



#' Subset PHMM
#' Modified from Cameron Nugent's COIL package https://github.com/CNuge/coil/blob/master/R/subset_PHMM.r
#' @param x A profile hidden Markov model (a "PHMM" object) generated by the aphid R package.
#' @param primers A character vector of length 2 contaning the forward and reverse primer sequences.
#' @param positions A numeric vector of length 2 containing the start and end position for model trimming. Used only if primers = NULL.
#' @param trimprimers logical indicating whether the primer-binding sites should be removed from the sequences in the returned model. Used only if positions = NULL.
#' @param min_score Minimum score for the viterbi alignment.
#' @param quiet Whether to print progress to console.
#'
#' @return
#' @export
#'
#' @examples
subset_model <- function(x, primers=NULL, positions=NULL, trimprimers=FALSE, min_score=10, quiet=FALSE){

  # Check inputs
  if (is.numeric(positions) & length(positions) == 2){
    start_bind <- sort(positions)[1]
    end_bind <- sort(positions)[2]
  } else if (is.null(positions) & (is.character(primers) & length(primers)== 2)){
    Fprimer <- primers[1]
    Rprimer <- primers[2]

    Fbind <-  get_binding_position(Fprimer, model=x, tryRC=TRUE, quiet=quiet, min_score = min_score)
    Rbind <-  get_binding_position(Rprimer, model=x, tryRC=TRUE, quiet=quiet, min_score = min_score)
    if(isTRUE(trimprimers)){
      start_bind <- Fbind$end + 1
      end_bind <- Rbind$start - 1
      if(!quiet){message("Primer binding positions will be trimmed from model")}
    } else if(isFALSE(trimprimers)){
      start_bind <- Fbind$start
      end_bind <- Rbind$end
    } else{
      stop("trimprimers must be either TRUE or FALSE")
    }

  } else if (is.null(positions) & is.null(primers)){
    stop("A value must be provided to either positions or primers")
  } else if(is.null(primers) & (!is.numeric(positions) | !length(positions) == 2) ){
    stop("positions must be a numeric vector of length 2 containing the start and end position for model trimming")
  } else if(is.null(positions) & (!is.character(primers) | !length(primers) == 2) ){
    stop("primers must be a character vector of length 2, with the forward primer as the first element, and reverse primer as the second element")
  } else if(!is.null(primers) & !is.null(positions)){
    stop("Only one of primers or positions can be entered at a time")
  } else {
    stop("A value must be provided to either positions or primers")
  }
  if(is.na(start_bind) & is.numeric(end_bind)){
    stop(paste0("Could not find binding position for Forward primer: ", Fprimer))
  } else if(is.na(end_bind) & is.numeric(start_bind)){
    stop(paste0("Could not find binding position for Reverse primer: ", Rprimer))
  } else if(is.na(start_bind) & is.na(end_bind)){
    stop(paste0("Could not find binding position for Forward primer: ", Fprimer, " and Reverse primer: ", Rprimer))
  }
  # Check that positions are right
  if(end_bind < start_bind){
    stop("Index error, end position must match or exceed the start position.")
  }
  if(end_bind > x$size){
    stop(paste0("Index error, end position out of bounds. Input PHMM only has a length of: ", x$size))
  }

  #get the map positions for subsetting the alignment-based fields
  map_start <- x$map[[start_bind]]
  map_end <- x$map[[end_bind]]

  #subset to get the new inserts
  new_inserts <- x$inserts[map_start:map_end]

  #subset to get the new insert lengths
  new_insert_lengths <- x$insertlengths[map_start:map_end]
  #reindex the insert length labels
  if(!is.null(new_insert_lengths)){
    names(new_insert_lengths) <- 0:(length(new_insert_lengths)-1)
  }
  #subset the mask
  new_mask <- x$mask[map_start:map_end]
  #subset the map
  new_map <- x$map[start_bind:end_bind]

  #set the first remaining map entry equal to one,
  #make the necessary relative adjustment to all other positions
  new_map <- new_map-(new_map[1]-1)
  #reindex the labels
  if(!is.null(new_map)){
    names(new_map) <- 1:(length(new_map))
  }

  #Subset the A and E matrixies, reindex the positions
  #A is 0 indexed, E is 1 indexed
  new_A <- x$A[,start_bind:end_bind]
  colnames(new_A) <- 0:(length(colnames(new_A))-1)
  new_E <- x$E[,start_bind:end_bind]
  colnames(new_E) <- 1:length(colnames(new_A))

  #Create the new PHMM object
  newPHMM = structure(list(name = x$name,
                           description = x$description,
                           size = (end_bind - start_bind)+1,
                           alphabet = x$alphabet,
                           A = new_A,
                           E = new_E,
                           qa = x$qa,
                           qe = x$qe,
                           inserts = new_inserts,
                           insertlengths = new_insert_lengths,
                           map = new_map,
                           date = x$date,
                           nseq = x$nseq,
                           weights = x$weights,
                           reference = x$reference,
                           mask = new_mask
  ), class =  "PHMM")
  message(paste0("PHMM model subset to ", newPHMM$size, " base pairs"))
  return(newPHMM)
}

# Alignment entropy -------------------------------------------------------

#' Alignment entropy
#'
#' @param x A DNAbin or AAbin object
#' @param mask_gaps The threshold of gaps allowed before a position in the alignment is masked
#' @param count_gaps Whether gaps should be counted within entropy calculations. Default is FALSE.
#' @param method 	the method employed by `entropy::entropy` to estimate alignment entropy. Accepts:
#' "ML" : maximum likelihood
#' "MM" : bias-corrected maximum likelihood,
#' "Jeffreys" : Dirichlet with a=1/2
#' "Laplace" : Dirichlet with a=1
#' "SG" : Dirichlet with a=a=1/length(y)
#' "minimax" : Dirichlet with a=sqrt(sum(y))/length(y)
#' "CS" : ChaoShen
#' "NSB": Nemenman, Shafee and Biale (2002)
#' "shrink" : Shrinkage estimator
#' See the help page of `entropy::entropy` for more information
#' @param unit the unit in which entropy is measured. The default is "nats" (natural units). For computing entropy in "bits" set unit="log2".
#' @param return_extra Whether to return a dataframe including extra columns including individual base counts, gap proportions and number of bases at each position
#'
#' @return
#' @export
#' @import purrr
#' @importFrom entropy entropy
#'
#'
#'
alignment_entropy <- function (x, mask_gaps=0.2, count_gaps = FALSE, method="ML", unit="log", return_extra=FALSE) {
  if ((mask_gaps < 0) | (mask_gaps > 1)) {
    stop("mask_gaps should be a percentage (within the [0,1] range).")
  }

  if(class(x)=="DNAbin"){
    x <- as.character(x)
    x <- lapply(x, toupper)
    names <- c("A", "C", "G", "T", "-")
  } else if(class(x)=="AAbin"){
    x <- as.character(x)
    x <- lapply(x, toupper)
    names <- c("A", "C", "D", "E", "F",
               "G", "H", "I", "K", "L",
               "M", "N", "P", "Q", "R",
               "S", "T", "V", "W", "Y", "-")
  }

  if(count_gaps){
    counter <- names
  } else if(!count_gaps){
    counter <- names[!names=="-"]
  }

  #Create matrix
  suppressWarnings(MSA <- matrix(as.vector(unlist(x)), ncol = length(x[[1]]), byrow = TRUE))

  n_seq <- length(MSA[,1])
  n_pos <- length(MSA[1,])

  #Summarise each position in alignment
  i=1
  df <- data.frame(matrix(ncol = length(names), nrow = n_pos))
  colnames(df) <- names

  for(i in 1:n_pos){
    df[i,names] <- sapply(names, function(x){ sum(MSA[, i] == x)})
  }
  ent <- df %>%
    dplyr::mutate(pos = rownames(.),
                  bases = n_seq - `-`,
                  prop_gaps = `-` / n_seq,
                  together = pmap(unname(.[counter]), c)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(ent = entropy::entropy(together, method=method, unit=unit)) %>%
    dplyr::mutate(ent = dplyr::case_when(
      prop_gaps > mask_gaps ~ as.numeric(NA),
      prop_gaps <= mask_gaps ~ ent
    )) %>%
    dplyr::select(-together)

  if(return_extra){
    out <- ent
  } else if(!return_extra){
    out <- ent$ent
    names(out) <- ent$pos
  }
  message("Masked ", sum(is.na(ent)), " alignment positions with over ", (mask_gaps*100), "% gaps")
  return(out)
}


# Summarise fasta ---------------------------------------------------------

#' summarise_fasta
#'
#' @param x The location of a fasta file or gzipped fasta file.
#' @param label optional, Add an extra column with a label
#' @param origin optional, a table with sequence id numbers and their database origins
#'
#' @return
#' @export
#' @import stringr
#' @import dplyr
#' @importFrom Biostrings fasta.index
#' @importFrom stats quantile
#'
#'
summarise_fasta <- function(x, label=NULL, origin=NULL) {
  if(is.null(origin)){
    out <- Biostrings::fasta.index(x) %>%
      dplyr::mutate(taxid = desc %>%
                      stringr::str_remove(pattern="(;)(.*?)(?=$)")  %>%
                      stringr::str_remove(pattern="(^)(.*?)(?<=\\|)")) %>%
      dplyr::summarise(nseqs = n(),
                       nspecies=n_distinct(taxid),
                       mean_length = mean(seqlength),
                       q0 = stats::quantile(seqlength, probs=0),
                       q25 = stats::quantile(seqlength, probs=0.25),
                       q50 = stats::quantile(seqlength, probs=0.5), # Q50 is median
                       q75 = stats::quantile(seqlength, probs=0.75),
                       q100 = stats::quantile(seqlength, probs=1)
      )

  } else if(is.data.frame(origin) | is_tibble(origin)){

    if(any(duplicated(origin$seqid))){
      stop("Origin table has duplicated seqids")
    }
    out <- Biostrings::fasta.index(x) %>%
      dplyr::mutate(taxid = desc %>%
                      stringr::str_remove(pattern="(;)(.*?)(?=$)")  %>%
                      stringr::str_remove(pattern="(^)(.*?)(?<=\\|)")) %>%
      dplyr::mutate(seqid = desc %>%
                      stringr::str_remove(pattern="(\\|)(.*?)(?=$)"))  %>%
      dplyr::left_join(origin, by="seqid") %>%
      dplyr::group_by(taxid) %>%
      dplyr::mutate(origin = paste(sort(unique(origin)), collapse="/"))%>% #Combine origins to make sure they arent duplicated
      dplyr::ungroup() %>%
      dplyr::group_by(origin) %>%
      dplyr::summarise(nseqs = n(),
                       nspecies=n_distinct(taxid),
                       mean_length = mean(seqlength),
                       q0 = stats::quantile(seqlength, probs=0),
                       q25 = stats::quantile(seqlength, probs=0.25),
                       q50 = stats::quantile(seqlength, probs=0.5), # Q50 is median
                       q75 = stats::quantile(seqlength, probs=0.75),
                       q100 = stats::quantile(seqlength, probs=1)
      )
  }
  if(is.character(label)) {
    out <- out %>%
      dplyr::mutate(label  = label)
  }
  return(out)
}


#' Download NCBI taxdump
#'
#' @param dest_dir A directory to save the downloaded ncbi taxdump files. If empty a new folder in the working directory will be created
#' @param include_synonyms Whether synonyms should be included in the database
#' @param force Whether already downloaded files should be overwrittenget
#'
#' @return
#' @export
#' @importFrom readr read_tsv
#' @importFrom utils untar
#' @importFrom utils download.file
#'
#'
get_ncbi_taxonomy <- function(dest_dir, include_synonyms = TRUE, force=FALSE) {
  if(missing(dest_dir)){
    dest_dir <- getwd()
  }
  tmp <- paste0(dest_dir,"/ncbi_taxdump")
  if (!dir.exists(tmp)) {
    dir.create(tmp) # Create first directory
  }
  if (force == TRUE | !file.exists(paste0(tmp, "/rankedlineage.dmp"))) {
    message("Downloading NCBI taxonomy database")
    fn <- "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
    utils::download.file(fn, destfile = paste0(tmp, "/tmp.tar.gz"))
    message("Extracting data\n")
    test <- utils::untar(tarfile = paste0(tmp, "/tmp.tar.gz"), exdir = tmp)
    if (!identical(test, 0L)) {
      stop(cat(test))
    }
    file.remove(paste0(tmp, "/tmp.tar.gz"))
  }
  message("Building NCBI taxonomy data frame\n")

  lin <- readr::read_tsv(paste0(tmp, "/rankedlineage.dmp"),
                         col_names = c("tax_id", "tax_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"),
                         col_types = ("i-c-c-c-c-c-c-c-c-c-")
  )

  # Remove synonyms
  if (!include_synonyms) {
    x <- scan(
      file = paste0(tmp, "/names.dmp"), what = "", sep = "\n",
      quiet = TRUE
    )
    syn <- x[grepl("synonym", x)]
    syn <- strsplit(syn, split = "\t")
    syn <- sapply(syn, function(s) s[c(1, 3)])
    syn <- as.data.frame(t(syn), stringsAsFactors = FALSE)
    syn[[1]] <- as.integer(syn[[1]])
    colnames(syn) <- c("taxID", "name")
    lin <- lin %>%
      filter(!tax_id %in% syn[[1]])
  }

  message("Done\n")
  attr(lin,'type') <- 'ncbi'
  return(lin)
}

# get_ncbi_lineage -------------------------------------------------------------
#' Get lineage
#'
#' @param x A DNAbin or DNAStringSet object with names in format `accession|tax_id;Genus species`
#' @param db a database file generated using ` get_ncbi_taxonomy()`
#'
#' @return
#' @export
#' @import stringr
#' @import dplyr
#' @importFrom tidyr separate
#' @importFrom tidyr unite
#' @importFrom tibble as_tibble
#'
get_ncbi_lineage <- function(x, db){
  if(missing(x)){stop("x must be a DNAbin or DNAStringSet object")}
  if(!length(x) > 0) {stop("x is empty or not a DNAbin or DNAStringSet object")}
  if(missing(db)){ db <-  get_ncbi_taxonomy()}
  cat("Getting taxonomic lineage from taxids\n")
  na_taxids <- names(x)[stringr::str_extract(names(x), "(?<=\\|).+?(?=;)") == "NA"]
  if(length(na_taxids)> 0){
    warning(length(na_taxids), " sequence/s have no matching tax_id in db")
  }
  lineage <- names(x) %>%
    tibble::as_tibble(column_name = "value") %>%
    tidyr::separate(col=value, into=c("acc", "species"), sep=";", extra="merge")%>%
    tidyr::separate(col=species, into=c("genus", "species"), sep=" |_", extra="merge")%>%
    tidyr::separate(col=acc, into=c("acc", "tax_id"), sep="\\|", extra="merge")%>%
    dplyr::mutate(tax_id = suppressWarnings(as.numeric(tax_id))) %>%
    dplyr::left_join (db %>% dplyr::select(-species, -genus), by = "tax_id")  %>%
    tidyr::unite(col = Acc, c("acc", "tax_id"), sep = "|")
  return(lineage)
}

# NCBI synonyms -----------------------------------------------------------
#' Get NCBI synonyms
#'
#' @param dir A directory containing the NCBI taxonomy that was downloaded using get_ncbi_taxonomy()
#' @param recurse Whether to recurse when searching for dir if dir is NULL.
#' If TRUE recurse fully, if a positive number the number of levels to recurse.
#' @param quiet Whether progress should be printed to the console.
#'
#' @return
#' @export
#' @import dplyr
#' @importFrom fs dir_ls
#'
#'
get_ncbi_synonyms <- function(dir=NULL, recurse=TRUE, quiet=FALSE) {
  if (is.null(dir)){
    if(!quiet){message("Dir is NULL, searching subdirectories for ncbi_taxdump")}
    dir <- fs::dir_ls(path = getwd(), type = "directory", glob = "*ncbi_taxdump", recurse = recurse)
    if(length(dir) >0){
      if(!quiet){message("ncbi_taxdump found at ", dir)}
    } else{
      stop("ERROR: provide a directory containing ncbi taxonomy")
    }
  }
  if(!quiet){message("Building synonyms data frame\n")}
  file <- normalizePath(paste0(dir, "/names.dmp"))

  #Read in
  splitLines <- do.call(rbind, strsplit(readLines(file),"\t\\|\t?"))
  splitLines<-splitLines[,-3]
  colnames(splitLines) <- c("tax_id","tax_name", "class")
  parsed <- as.data.frame(splitLines)

  out <- parsed %>%
    dplyr::group_by(tax_id) %>%
    dplyr::filter(any(class=="synonym")) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!class=="scientific name")  %>%
    dplyr::select(-class, synonym = tax_name) %>%
    dplyr::left_join(parsed %>%
                  dplyr::filter(class=="scientific name") %>%
                  dplyr::select(tax_id, tax_name), by="tax_id")

  return(out)
}

# resolve_synonyms_ncbi ---------------------------------------------------
#' resolve_synonyms_ncbi
#'
#' @param x A DNAbin or DNAStringSet Object
#' @param dir A directory containing the NCBI taxonomy that was downloaded using get_ncbi_taxonomy()
#' @param quiet Whether progress should be printed to console
#'
#' @return
#' @export
#' @import stringr
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom tidyr separate
#' @importFrom tidyr unite
#' @importFrom ape as.DNAbin
#' @importFrom methods is
#'
resolve_synonyms_ncbi <- function(x, dir=NULL, quiet = FALSE) {
  time <- Sys.time() # get time
  if (quiet == TRUE) {
    verbose <- FALSE
  } else {
    (verbose <- TRUE)
  }

  # Check type of input
  if (methods::is(x, "DNAbin")) {
    message("Input is DNAbin")
    is.seq <- TRUE
  } else  if (methods::is(x, "DNAStringSet")| methods::is(x, "DNAString")) {
    message("Input is DNAStringSet, converting to DNAbin")
    x <- ape::as.DNAbin(x)
    is.seq <- TRUE
  } else  if (methods::is(x, "character")) {
    message("Input is character vector, resolving Genus species binomials")
    is.seq <- FALSE
  }

  syns <- get_ncbi_synonyms(dir=dir) #...=...

  # if input has sequences, get names
  if (is.seq) {
    query <- names(x) %>%
      stringr::str_split_fixed(";", n = 2) %>%
      as.data.frame() %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(query = V2) %>%
      dplyr::mutate(query = stringr::str_replace_all(query, "_", " "))
  } else if (!is.seq) {
    query <- data.frame(query= x %>% stringr::str_replace_all(query, "_", " "))
    query$tax_id <- "NA" # Add dummy columns
    query$acc <- "NA"
  }

  #Get synonyms to replace
  to_replace <- query %>%
    dplyr::select(-tax_id) %>%
    dplyr::filter(query %in% syns$synonym) %>%
    dplyr::left_join(syns %>% dplyr::rename(query = synonym), by="query") %>%
    dplyr::select(-query) %>%
    dplyr::filter(!duplicated(acc)) #where are these coming from?

  if(nrow(to_replace) > 0){
  out <- query %>%
    dplyr::rename(tax_name = query) %>%
    dplyr::rows_update(to_replace, by="acc") %>%
    tidyr::unite(col = V1, c("acc", "tax_id"), sep = "|") %>%
    dplyr::mutate(tax_name = tax_name %>% stringr::str_replace_all(" ", "_"))

  if(is.seq) {
      names(x) <- paste(out$V1, out$tax_name, sep = ";")
    } else if (!is.seq) {
      x  <- paste(out$V1, out$tax_name, sep = ";")
    }
  time <- Sys.time() - time
  if (!quiet) {message(paste0("resolved ", nrow(to_replace), " synonyms in ", format(time, digits = 2)))}
  } else (message("No synonyms detected"))
  return(x)
}

# ncbi_taxid --------------------------------------------------------------
#' Get ncbi taxid's for a taxon name
#'
#' @param x A DNAbin or DNAStringSet object with names in format `accession|tax_id;Genus species`
#' @param db a database file generated using ` get_ncbi_taxonomy()`
#'
#' @return
#' @export
#' @importFrom magrittr set_colnames
#' @import dplyr
#'
#'
ncbi_taxid <- function(x, db=NULL) {

  if (is.null(db)) { db <- get_ncbi_taxonomy()}
  out <-  as.data.frame(x) %>%
    magrittr::set_colnames("tax_name") %>%
    dplyr::left_join (db, by = "tax_name") %>%
    dplyr::pull(tax_id)

  return(out)
}


# gid_to_acc --------------------------------------------------------------
#' Convert NCBI gene ids to accession numbers
#'
#' @param ids A character or numeric vector of NCBI gids
#' @param database The origin database of the gids, default 'nuccore'
#' @param chunk_size The size of the chunked searches to conduct.
#' Warning, chunk sizes over 300 can be too big for the NCBI servers.
#' @param multithread Whether multithreading should be used
#' @param progress Whether a progress bar is displayed.
#' @param quiet Whether progress should be printed to the console.
#'
#' @return
#' @export
#' @import future
#' @import furrr
#' @import purrr
#' @importFrom rentrez entrez_post
#' @importFrom rentrez entrez_fetch
#'
#'
gid_to_acc <- function(ids, database="nuccore", chunk_size=300, multithread=TRUE, progress=FALSE, quiet=FALSE){
  if(!class(ids) %in% c("character", "numeric")){
    stop("input ids must be a character or numeric vector of NCBI gids")
  }
  time <- Sys.time() # get time

  # split search into chunks
  chunks <- split(ids, ceiling(seq_along(ids)/chunk_size))
  if(!quiet){message(paste0("Converting ", length(ids), " NCBI gids to accession numbers in ", length(chunks), " chunks"))}

  # setup multithreading
  cores <- setup_multithread(multithread)

  #Main function
  out <- furrr::future_map(chunks, function(x){
    upload <- rentrez::entrez_post(db=database, id=x)
    dl <- rentrez::entrez_fetch(db = database, web_history = upload, rettype = "acc", retmax = chunk_size)
    acc <- readLines(textConnection(dl))
    acc <- acc[!acc==""]
    return(acc)
  }, .progress = progress) %>%
    unlist(use.names=FALSE)

  #finished
  time <- Sys.time() - time
  if (!quiet) {message(paste0("Sucessfully converted ", length(out), " of ", length(ids), " NCBI gids to accession numbers in ", format(time, digits = 2)))}
  return(out)
}


# Get_ott_taxonomy --------------------------------------------------------

#' Download open tree of life taxonomy
#'
#' @param url a URL to download from, if left blank the latest version will be downloaded
#' @param dest_dir A directory to save the zipped taxonomy database to, if left blank the working directory is selected
#' @param force Whether existing files should be overwritten
#'
#' @return
#' @export
#' @import stringr
#' @importFrom httr GET
#' @importFrom httr write_disk
#' @importFrom utils untar
#' @importFrom rvest html_nodes
#' @importFrom rvest html_attr
#' @importFrom xml2 read_html
#'
download_ott_taxonomy <- function(url, dest_dir, force=FALSE) {

  if (missing(url)) {
    # find the latest version of taxonomy
    download_page <- xml2::read_html("https://tree.opentreeoflife.org/about/taxonomy-version/")
    link_hrefs <- download_page %>%
      rvest::html_nodes("a") %>%
      rvest::html_attr("href")
    url <- grep("http://files.opentreeoflife.org/ott/.*tgz$", link_hrefs, perl = TRUE) %>%
      link_hrefs[.] %>% .[1]
  }

  version <- basename(url) %>%
    stringr::str_remove(".tgz")

  if(missing(dest_dir)){
    message("dest_dir is missing, downloading into /", version)
    dest_dir <- version
  }

  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir) # Create first directory
  }

  # Check if files exist already in dest_dir
  expected_files <- c(
    paste0(dest_dir,"/", version,".tgz" ),
    paste0(dest_dir,"/taxonomy.tsv" ),
    paste0(dest_dir,"/synonyms.tsv" )
  )

  if (any(file.exists(expected_files)) & force == FALSE) {
    message(paste0("Skipped as files already exist in dest_dir, to overwrite set force to TRUE"))
    return(NULL)
  } else  if (any(file.exists(expected_files)) && force == TRUE) {
    unlink(expected_files, recursive = TRUE) # Remove existing version
  }

  dest_file <- file.path(dest_dir, basename(url))
  httr::GET(url, httr::write_disk(dest_file, overwrite=TRUE))

  #unzip file
  utils::untar(dest_file, exdir = ".")
  #Remove download
  file.remove(dest_file)
  message(paste0("Downloaded taxonomy to: ",stringr::str_remove(dest_file,".tgz" ), " \n"))
  return(stringr::str_remove(dest_file,".tgz" ))
}


# get_ott_taxonomy ------------------------------------------------------

#' get_ott_taxonomy
#'
#' @param dir a directory containing the open tree of life taxonomy files obtained from the  `download_ott_taxonomy` function
#' @param quiet Whether progress should be printed to console
#' @param filter_unplaced Whether to filter 'bad' entries. These include
#' incertae_sedis
#' major_rank_conflict
#' unplaced
#' environmental
#' inconsistent
#' extinct
#' hidden
#' hybrid
#' not_otu
#' viral
#' barren
#' See: https://github.com/OpenTreeOfLife/reference-taxonomy/blob/master/doc/taxon-flags.md for more info
#'
#' @return
#' @export
#' @import dplyr
#' @importFrom vroom vroom
#' @importFrom data.table data.table
#' @importFrom data.table tstrsplit
#'
get_ott_taxonomy <- function(dir=NULL, quiet=FALSE, filter_unplaced=TRUE) {
  if (is.null(dir)){
    input <- NA
    while(!isTRUE(input == "1") && !isTRUE(input == "2")) {
      input <- readline(prompt="No directory provided, type '1' if you want to download the ott taxonomy \n")
      if(input == "1") {
        dir <- download_ott_taxonomy()
      } else {stop("Stopped function")}
    }
  }
  if(!quiet){message("Building data frame\n")}
  file <- normalizePath(paste0(dir, "/taxonomy.tsv"))

  if (filter_unplaced==TRUE){
    remap <- vroom::vroom(file, delim="\t|\t") %>%
      dplyr::filter(!grepl("incertae_sedis|incertae_sedis$|major_rank_conflict|unplaced|environmental|inconsistent|extinct|hidden|hybrid|not_otu|viral|barren", flags)) %>%
      dplyr::mutate(rank = stringr::str_replace(rank, pattern="no rank - terminal", replacement="terminal")) %>%
      dplyr::select(uid, parent_uid, name, rank, sourceinfo, flags) %>%
      dplyr::rename(tax_id = uid, parent_taxid = parent_uid, tax_name = name)
  }else {
    remap <- vroom::vroom(file, delim="\t|\t") %>%
      dplyr::select(uid, parent_uid, name, rank, sourceinfo, flags) %>%
      dplyr::rename(tax_id = uid, parent_taxid = parent_uid, tax_name = name)
  }
  d.dt <- data.table::data.table(remap, key="tax_id")
  db <- d.dt[, list(sourceinfo = unlist(strsplit(sourceinfo, ",")), tax_name, parent_taxid, rank, flags), by=tax_id
  ][, c("source", "id") := data.table::tstrsplit(sourceinfo, ":", fixed=TRUE)
  ][,c('sourceinfo') :=  .(NULL)]
  attr(db,'type') <- 'OTT'
  return(db)
}

# map_to_ott  ------------------------------------------------------------

#' Map taxa to open tree of life
#'
#' @param x a DNAbin
#' @param db an OTT taxonomic database
#' @param from the existing taxonomic ID format
#' @param resolve_synonyms Whether to resolve synonyms
#' @param dir A directory containing the OTT taxonomy, required if resolve synonyms is true
#' @param filter_unplaced Whether to filter 'bad' entries. These include
#' incertae_sedis
#' major_rank_conflict
#' unplaced
#' environmental
#' inconsistent
#' extinct
#' hidden
#' hybrid
#' not_otu
#' viral
#' barren
#' See: https://github.com/OpenTreeOfLife/reference-taxonomy/blob/master/doc/taxon-flags.md for more info
#' @param remove_na Whether taxa that could not be mapped to the open tree of life should be removed
#' @param quiet Whether progress should be printed to console
#'
#' @return
#' @export
#' @import dplyr
#' @import stringr
#' @importFrom tidyr separate
#' @importFrom tibble as_tibble
#' @importFrom ape as.DNAbin
#' @importFrom methods is
#'
map_to_ott <- function(x, db, from="ncbi", resolve_synonyms=TRUE, dir=NULL, filter_unplaced=TRUE, remove_na = FALSE, quiet=FALSE){
  time <- Sys.time() # get time

  #input checks
  if(resolve_synonyms && is.null(dir)){
    stop("If resolve_synonmys is true, a directory containing the OTT taxonomy must be provided")
  }
  if(!attr(db,'type')=="OTT") {stop("Error: requires OTT db, generate one with get_ott_taxonomy")}
  if(missing(db) & is.null(dir)) {stop("Error: requires OTT db, generate one with get_ott_taxonomy")}
  if(missing(db) & file.exists(paste0(dir,"/taxonomy.tsv"))){
    message("No db provided but dir provided, creating new db")
    db <- get_ott_taxonomy(dir=dir, filter_unplaced = filter_unplaced)
  }
  if (!from %in% unique(db$source)){ stop("Error: 'from' is not in db")}

  #Check input format
  if (methods::is(x, "DNAbin")) {
    message("Input is DNAbin")
    tax <- names(x) %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else if (methods::is(x, "DNAStringSet")) {
    message("Input is DNAStringSet")
    tax <- ape::as.DNAbin(x) %>%
      names() %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  }else  if (methods::is(x, "character") && (str_detect(x, "\\|") & str_detect(x, ";"))) {
    message("Detected | and ; delimiters, assuming 'Accession|taxid;Genus Species' format")
    tax <- x %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else  if (methods::is(x, "character") && !(stringr::str_detect(x, "\\|") && stringr::str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of species names")
    tax <- data.frame(acc = as.character(NA), id=as.character(NA), tax_name = x, stringsAsFactors = FALSE)
  }else (stop("x must be DNA bin or character vector"))

  #Read in synonyms DB
  #TODO: make synonyms df an attribute of the DB object
  if(resolve_synonyms == TRUE){ syn <- parse_ott_synonyms(dir=dir)}

  #Get tax
  tax <- tax %>%
    dplyr::left_join (db %>%   # First map by id
                        dplyr::filter(source==!!from) %>%
                        dplyr::select(-source) %>%
                        dplyr::rename(tax_name.x = tax_name), by = "id") %>%
    dplyr::left_join (db %>%   # Next map by tax_name
                        dplyr::select(-id, -source ) %>%
                        dplyr::filter(!duplicated(tax_name)), by = "tax_name") # then map by name

  #Resolve synonyms
  if(resolve_synonyms == TRUE){
    tax <- tax %>%
      dplyr::left_join(syn %>%
                         dplyr::rename(tax_name.z = tax_name, tax_name = synonym, tax_id.z = tax_id) %>%
                         dplyr::filter(!duplicated(tax_name)),
                       by = "tax_name") %>%
      dplyr::mutate(tax_id = dplyr::case_when(
        !is.na(tax_id.z) ~ tax_id.z, #If synonym was found, use it
        !is.na(tax_id.x) & is.na(tax_id.z) ~ tax_id.x, #If no synonym was found, but an ID match was found use it
        is.na(tax_id.x) & is.na(tax_id.z) & !is.na(tax_id.y)  ~ tax_id.y #If no synonym and no ID match found, but a name match was, use it
      ),
      tax_name = dplyr::case_when(
        !is.na(tax_name.z) ~ tax_name.z, #If synonym was found, use it
        !is.na(tax_name.x) & is.na(tax_name.z) ~ tax_name.x, #If no synonym was found, but an ID match was found use it
        is.na(tax_name.x) & is.na(tax_name.z) ~ tax_name  #If no synonym was found, and no ID match, retain current name
      ))

    if (filter_unplaced == TRUE){ #ensure resolving synonyms didnt introduce bads
      bads <- db %>%
        dplyr::filter(grepl("incertae_sedis,|incertae_sedis$|major_rank_conflict|unplaced|environmental|inconsistent|extinct|hidden|hybrid|not_otu|viral|barren", flags))
      tax <- tax %>%
        dplyr::mutate(tax_id = dplyr::case_when( #Ensure no
          tax_name %in% bads$name ~  as.numeric(NA),
          !tax_name %in% bads$name ~ tax_id
        )) %>%
        dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))
    } else if (filter_unplaced == FALSE){
      tax <- tax %>%
        dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))
    }
  } else if( resolve_synonyms == FALSE){
    tax <- tax %>%
      dplyr::mutate(tax_id = dplyr::case_when(
        !is.na(tax_id.x) ~ tax_id.x,
        is.na(tax_id.x) & !is.na(tax_id.y) ~ tax_id.y
      ),
      tax_name = dplyr::case_when(
        is.na(tax_name.x) ~ tax_name,
        !is.na(tax_name.x) ~ tax_name.x
      )) %>%
      dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))
  }

  #Replace names
  if (methods::is(x, "DNAbin") | methods::is(x, "DNAStringSet")) {
    names(x) <- tax$name
  } else  if (methods::is(x, "character")) {
    x <- tax$name
  }
  time <- Sys.time() - time
  if (!quiet) (message(paste0("Translated ",  length(x)," tax_ids from ",from, " to Open tree of life in ", format(time, digits = 2))))

  # Filter NA's
  if (remove_na ==TRUE){
    remove <- tax %>%
      dplyr::filter(is.na(tax_id)) %>%
      dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))

    if(methods::is(x, "DNAbin") | methods::is(x, "DNAStringSet")){
      x[names(x) %in% remove$name] <- NULL
    }else if (methods::is(x, "character")){
      x[x %in% remove$name] <- NULL
    }
    if(!quiet){message(paste0("Removed ", nrow(remove), " sequences that could not be mapped to OTT\n"))}
  }
  return(x)
}

# Parse Synonyms  -----------------------------------------------------------

#' parse the open tree of life synonyms file
#'
#' @param dir a directory containing the open tree of life taxonomy files obtained from the  `download_ott_taxonomy` function
#' @param quiet Whether progress should be printed to console
#'
#' @return
#' @export
#' @import stringr
#' @import dplyr
#' @importFrom vroom vroom
#'
#'
parse_ott_synonyms <- function(dir=NULL, quiet=FALSE) {
  if (is.null(dir)){
    stop("ERROR: provide a directory containing ott taxonomy")
  }
  if(!quiet){message("Building synonyms data frame\n")}
  file <- normalizePath(paste0(dir, "/synonyms.tsv"))
  out <- vroom::vroom(file, delim="\t|\t" )%>%
    dplyr::mutate(tax_name = uniqname %>%
                    stringr::str_remove(pattern="^.*for ")%>%
                    stringr::str_remove(pattern="\\).*$")%>%
                    stringr::str_remove(pattern="\\(.*$")
    ) %>%
    dplyr::rename(tax_id = uid, synonym = name) %>%
    dplyr::select(tax_id, tax_name, synonym)
  return(out)
}

# OTT Recursion -----------------------------------------------------------


#' Recursively get lineage from OTT taxid
#' This function derives the full lineage of a taxon ID number from a given taxonomy database
#' @param x A DNAbin, DNAStringSet, taxonomy headers formnatted Acc|taxid;taxonomy, or vector of taxids
#' @param db an OTT taxonomic database
#' @param ranks the taxonomic ranks to filter to. Default is "kingdom", "phylum", "class", "order", "family", "genus", "species"
#' To get strain level ranks, add "terminal" to ranks
#' @param output Options are "tax_name" to return taxonomic names for each rank in the lineage or "tax_id" to return taxonomic ID's.
#' @param multithread Whether multithreading should be used, if TRUE the number of cores will be automatically detected (Maximum available cores - 1), or provide a numeric vector to manually set the number of cores to use. Default is FALSE (single thread)
#' This argument may alternatively be a 'cluster' object, in which case it is the user's responsibility to close the socket connection at the conclusion of the operation,
#' for example by running parallel::stopCluster(cores).
#'
#' @return
#' @export
#' @import parallel
#' @import dplyr
#' @import stringr
#' @importFrom tidyr separate
#' @importFrom tidyr unite
#' @importFrom tidyr pivot_wider
#' @importFrom tibble as_tibble
#' @importFrom methods is
#'
get_ott_lineage <- function(x, db, output="tax_name", ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), multithread = FALSE){
  #Check input format
  if (methods::is(x, "DNAbin")) {
    message("Input is DNAbin")
    lineage <- names(x) %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else if (methods::is(x, "DNAStringSet")) {
    message("Input is DNAStringSet")
    lineage <- as.DNAbin(x) %>%
      names() %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  }else  if (methods::is(x, "character") && (stringr::str_detect(x, "\\|") & stringr::str_detect(x, ";"))) {
    message("Detected | and ; delimiters, assuming 'Accession|taxid;Genus Species' format")
    lineage <- x %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else  if (methods::is(x, "character") && !(stringr::str_detect(x, "\\|") && stringr::str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of tax_id's")
    lineage <- data.frame(acc = as.character(NA), tax_name=as.character(NA), tax_id = x, stringsAsFactors = FALSE)
  }else (stop("x must be DNA bin or character vector"))

  # setup multi threading
  cores <- setup_para(multithread, quiet)

  # Check for duplicated accessions
  if(any(duplicated(lineage$acc))){stop("Duplicated sequence accessions found")}

  #check db
  if(missing(db) | !attr(db,'type')=="OTT") {stop("Error: requires OTT db, generate one with get_ott_taxonomy")}

  tax_ids <- as.numeric(lineage$tax_id)
  db$rank <- as.character(db$rank)
  db$tax_name <- as.character(db$tax_name) # avoid stringsasfactor issues

  #dereplicate to uniques
  uh <- unique(paste(tax_ids))
  pointers <- seq_along(uh)
  names(pointers) <- uh
  pointers <- unname(pointers[paste(tax_ids)])
  tax_ids <- tax_ids[!duplicated(pointers)]

  #Recursive function
  gl1 <- function(tax_id, db, ranks){
    if(is.na(tax_id)) return(NA)
    stopifnot(length(tax_id) == 1 & mode(tax_id) == "numeric")
    res <- character(100)
    resids <- integer(100)
    resnames <- character(100)
    counter <- 1
    index <- match(tax_id, db$tax_id)
    if(is.na(index)){
      # warning(paste("Taxon ID", tax_id, "not found in database\n"))
      return(data.frame(rank= NA,
                        tax_name=NA,
                        tax_id=NA,
                        stringsAsFactors = FALSE))
    }
    repeat{
      if(is.na(index)) break
      # if(length(index) > 1) cat(index, "\n")
      res[counter] <- db$tax_name[index]
      resids[counter] <- db$tax_id[index]
      resnames[counter] <- db$rank[index]
      index <- db$parent_tax_index[index]
      counter <- counter + 1
    }
    #get position of ranks
    pos <- match(ranks, resnames)
    resnames <- resnames[pos]
    res <- res[pos]
    resids <- resids[pos]
    out <- data.frame(rank= ranks,
                      tax_name=res,
                      tax_id=resids,
                      stringsAsFactors = FALSE)
    return(out)
  }
  db$parent_tax_index <- match(db$parent_taxid, db$tax_id)
  ## multithreading
  if(inherits(cores, "cluster")){
    res <- parallel::parLapply(cores, tax_ids, gl1, db, ranks)
  }else if(cores == 1){
    res <- lapply(tax_ids, gl1, db, ranks)
  }else{
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' argument")
    if(cores > 1){
      cl <- parallel::makeCluster(cores)
      res <- parallel::parLapply(cl, tax_ids, gl1, db, ranks)
      parallel::stopCluster(cl)
    }else{
      res <- lapply(tax_ids, gl1, db, ranks)
    }
  }
  res <- res[pointers] #re-replicate
  names(res) <- lineage$acc
  if(output =="tax_name"){
    out <- bind_rows(res, .id="id") %>%
      dplyr::select(-tax_id) %>%
      dplyr::group_by(id) %>%
      tidyr::pivot_wider(
        names_from = rank,
        values_from = tax_name) %>%
      dplyr::ungroup() %>%
      dplyr::bind_cols(lineage) %>%
      tidyr::unite(Acc, c(acc, tax_id), sep = "|") %>%
      dplyr::select(Acc, all_of(ranks), tax_name)
  } else if(output == "tax_id"){
    out <- bind_rows(res, .id="id") %>%
      dplyr::select(-tax_name) %>%
      dplyr::group_by(id) %>%
      tidyr::pivot_wider(
        names_from = rank,
        values_from = tax_id) %>%
      dplyr::ungroup() %>%
      dplyr::bind_cols(lineage) %>%
      tidyr::unite(Acc, c(acc, tax_id), sep = "|") %>%
      dplyr::select(Acc, all_of(ranks), tax_name)
  }
  return(out)
}


# Filter ott ---------------------------------------------------------

#' Filter unplaced taxonomic labels
#' @description
#' Filter flags indicating unplaced taxa in the taxonomic tree. These include:
#' incertae_sedis
#' major_rank_conflict
#' unplaced
#' environmental
#' inconsistent
#' extinct
#' hidden
#' hybrid
#' not_otu
#' viral
#' barren
#' See: https://github.com/OpenTreeOfLife/reference-taxonomy/blob/master/doc/taxon-flags.md for more info
#'
#' @param x A DNAbin, DNAStringSet, taxonomy headers formnatted Acc|taxid;taxonomy, or vector of taxids
#' @param db an OTT taxonomic database
#' @param quiet Whether progress should be printed to the console.
#'
#' @return
#' @export
#' @import stringr
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom tidyr separate
#' @importFrom ape as.DNAbin
#' @importFrom methods is
#'
filter_unplaced <- function(x, db, quiet=FALSE){
  #Check input format
  if (methods::is(x, "DNAbin")) {
    message("Input is DNAbin")
    tax <- names(x) %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else if (methods::is(x, "DNAStringSet")) {
    message("Input is DNAStringSet")
    tax <- ape::as.DNAbin(x) %>%
      names() %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  }else  if (methods::is(x, "character") && (str_detect(x, "\\|") & str_detect(x, ";"))) {
    message("Detected | and ; delimiters, assuming 'Accession|taxid;Genus Species' format")
    tax <- x %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else  if (methods::is(x, "character") && !(str_detect(x, "\\|") && str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of species names")
    tax <- data.frame(acc = as.character(NA), tax_id=as.character(NA), tax_name = x, stringsAsFactors = FALSE)
  }else (stop("x must be DNA bin or character vector"))

  #check db
  if(missing(db) | !attr(db,'type')=="OTT") {stop("Error: requires OTT db, generate one with get_ott_taxonomy")}
  bads <- db %>%
    dplyr::filter(grepl("incertae_sedis,|incertae_sedis$|major_rank_conflict|unplaced|environmental|inconsistent|extinct|hidden|hybrid|not_otu|viral|barren", flags))
  if(nrow(bads)==0){stop("Database is already filtered for bad flags, rerun get_ott_taxonomy() with filter_unplaced=FALSE")}
  tax <- tax %>%
    dplyr::mutate(tax_id = dplyr::case_when( #Ensure no
      tax_id %in% bads$tax_id ~  as.character(NA),
      !tax_id %in% bads$tax_id ~ tax_id
    )) %>%
    dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))

  # Filter NA's
  remove <- tax %>%
    dplyr::filter(is.na(tax_id)) %>%
    dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))

  if(methods::is(x, "DNAbin") | methods::is(x, "DNAStringSet")){
    names(x) <- tax$name
    x <- x[!names(x) %in% remove$name]

  }else if (methods::is(x, "character")){
    x <- tax$name
    x <- x[!x %in% remove$name]
  }
  if(!quiet){message(paste0("Removed ", nrow(remove), " unplaced flagged taxa\n"))}
  return(x)

}

#' Filter infraspecific taxonomic labels
#' @param x A DNAbin, DNAStringSet, taxonomy headers formnatted Acc|taxid;taxonomy, or vector of taxids
#' @param db an OTT taxonomic database
#' @param quiet Whether progress should be printed to the console.
#'
#' @return
#' @export
#' @import stringr
#' @importFrom tidyr separate
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom methods is
#'
filter_infraspecifc <- function(x, db, quiet=FALSE){
  #Check input format
  if (methods::is(x, "DNAbin")) {
    message("Input is DNAbin")
    tax <- names(x) %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else if (methods::is(x, "DNAStringSet")) {
    message("Input is DNAStringSet")
    tax <- ape::as.DNAbin(x) %>%
      names() %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  }else  if (methods::is(x, "character") && (str_detect(x, "\\|") & str_detect(x, ";"))) {
    message("Detected | and ; delimiters, assuming 'Accession|taxid;Genus Species' format")
    tax <- x %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else  if (methods::is(x, "character") && !(str_detect(x, "\\|") && str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of species names")
    tax <- data.frame(acc = as.character(NA), tax_id=as.character(NA), tax_name = x, stringsAsFactors = FALSE)
  }else (stop("x must be DNA bin or character vector"))

  #check db
  if(missing(db) | !attr(db,'type')=="OTT") {stop("Error: requires OTT db, generate one with get_ott_taxonomy")}
  bads <- db %>%
    dplyr::filter(grepl("infraspecific", flags))
  tax <- tax %>%
    dplyr::mutate(tax_id = dplyr::case_when( #Ensure no
      tax_id %in% bads$tax_id ~  as.character(NA),
      !tax_id %in% bads$tax_id ~ tax_id
    )) %>%
    dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))

  # Filter NA's
  remove <- tax %>%
    dplyr::filter(is.na(tax_id)) %>%
    dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))

  if(methods::is(x, "DNAbin") | methods::is(x, "DNAStringSet")){
    names(x) <- tax$name
    x <- x[!names(x) %in% remove$name]

  }else if (methods::is(x, "character")){
    x <- tax$name
    x <- x[!x %in% remove$name]
  }
  if(!quiet){message(paste0("Removed ", nrow(remove), " infraspecific taxa\n"))}

  # check for multiple taxids per name remaining
  checks <- tax %>%
    dplyr::filter(!is.na(tax_id)) %>%
    dplyr::group_by(tax_name) %>%
    dplyr::summarise(count = n_distinct(tax_id))

  if (any(checks$count >1 )) {warning("Multiple tax_ids remain for some tax_names, please check manually")}
  return(x)
}



# Reformat hierarchy ------------------------------------------------------

#' Reformat hierarchy
#'
#' @param x A DNAbin with names formatted Accession|taxid;taxonomy
#' @param db A database file generated by `get_ncbi_taxonomy` or `get_ott_lineage`. Required if get_lineage is TRUE
#' @param quiet Whether progress should be printed to console
#' @param ranks The taxonomic ranks to return in the output names
#' @param multithread Whether multithreading should be used, if TRUE the number of cores will be automatically detected (Maximum available cores - 1), or provide a numeric vector to manually set the number of cores to use. Default is FALSE (single thread)
#' This argument may alternatively be a 'cluster' object, in which case it is the user's responsibility to close the socket connection at the conclusion of the operation,
#' for example by running parallel::stopCluster(cores).
#'
#' @return
#' @export
#' @import stringr
#' @import dplyr
#' @importFrom tidyr separate
#' @importFrom tidyr unite
#' @importFrom magrittr set_colnames
#' @importFrom tibble as_tibble
#' @importFrom ape as.DNAbin
#' @importFrom methods is
#'
reformat_hierarchy <- function(x, db = NULL, quiet = FALSE,
                               ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), multithread = FALSE) {
  time <- Sys.time() # get time
  # Convert to DNAbin
  if (!methods::is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
  }
  if (missing(db)) {
    stop("A taxonomic database needs to be provided, generate one with get_ncbi_taxonomy or get_ott_taxonomy")
  }
  if(attr(db, "type") == "ncbi"){
    # Split current names
    seqnames <- names(x) %>%
      stringr::str_remove(";$") %>%
      stringr::str_split_fixed(";", n = Inf) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|")%>%
      dplyr::select(1:2, tail(names(.), 1)) %>%
      magrittr::set_colnames(c("acc", "tax_id", "species")) %>%
      dplyr::mutate(tax_id = suppressWarnings(as.numeric(tax_id)))

    # Get lineage from taxid
    lineage <- seqnames %>%
      dplyr::left_join (db %>% dplyr::select(-species), by = "tax_id")  %>%
      tidyr::unite(col = Acc, c(acc, tax_id), sep = "|") %>%
      tidyr::unite(col = names, c(!!ranks), sep = ";")
  } else if(attr(db, "type") == "OTT"){
    lineage <- get_ott_lineage(x, db=db, ranks=ranks, multithread=multithread) %>%
      tidyr::unite(col = names, c(!!ranks), sep = ";")
  }

  names(x) <- paste(lineage$Acc, lineage$names, "", sep = ";") %>%
  stringr::str_remove(";$")
  time <- Sys.time() - time
  if (!quiet) (message(paste0("Reformatted annotations for ", length(x), " sequences in ", format(time, digits = 2))))
  return(x)
}


# Reformat DADA2 Genus ------------------------------------------------------
#' Reformat DADA2 Genus
#' @param x A DNAbin or DNAStringSet object with names in format `accession|tax_id;Genus species`
#' @param quiet Whether progress should be printed to the console.
#' @param ranks The taxonomic ranks to include in output. default is kingdom;phylum;class;order;family;genus
#' @param db A database file generated by `get_ncbi_taxonomy` or `get_ott_lineage`. Required if get_lineage is TRUE
#'
#' @return
#' @export
#'
#'
reformat_dada2_gen <- function(x, db = NULL, ranks = c("kingdom", "phylum", "class", "order", "family", "genus"), quiet = FALSE) {
  reformat_hierarchy(x, db = db, quiet = quiet, ranks = ranks)
}


# Reformat DADA2 Species ----------------------------------------------------
#' Reformat DADA2 Species
#'
#' @param x A DNAbin with names formatted Accession|taxid;taxonomy
#' @param quiet Whether progress should be printed to console
#'
#' @return
#' @export
#' @import stringr
#' @import dplyr
#' @importFrom magrittr set_colnames
#' @importFrom tibble as_tibble
#' @importFrom ape as.DNAbin
#'
reformat_dada2_spp <- function(x, quiet = FALSE) {
  time <- Sys.time() # get time
  # Convert to DNAbin
  if (!methods::is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
  }
  # Split current names
  seqnames <- names(x) %>%
    stringr::str_replace(";$", "") %>%
    stringr::str_split_fixed(";", n = Inf) %>%
    tibble::as_tibble() %>%
    dplyr::select(1, tail(names(.), 1)) %>%
    magrittr::set_colnames(c("acc", "species")) %>%
    dplyr::mutate(species = stringr::str_replace(species, pattern="_", replacement=" "))

  names(x) <- paste(seqnames$acc, seqnames$species, sep = " ")
  time <- Sys.time() - time
  if (!quiet) (message(paste0("Reformatted annotations for ", length(x), " sequences in ", format(time, digits = 2))))
  return(x)
}


# train IDTAXA -----------------------------------------------------------

#' Train IDTAXA
#'
#' @param x A DNAbin
#' @param max_group_size The maximum size of any taxonomic group. This can be set to Inf (infinity) to allow for an unlimited number of sequences per group.
#' @param max_iterations The number of iterations to conduct training for in order to identify  any training sequences whose assigned classifications completely disagree with their predicted classification.
#' @param allow_group_removal  Whether sequences that are the last remaining representatives of an entire group in the training data can be removed.
#' This can occur if the entire group appears to be misplaced in the taxonomic tree.
#' @param orient Training sequences must all be in the same orientation. Set this to TRUE to reorient the sequences if you are unsure.
#' @param get_lineage Get full taxonomic lineage using reformat_hierarchy if not already present.
#' @param db A database file generated by `get_ncbi_taxonomy` or `get_ott_lineage`. Required if get_lineage is TRUE.
#' @param quiet Whether progress should be printed to console.
#'
#' @return
#' @export
#' @import stringr
#' @importFrom DECIPHER RemoveGaps
#' @importFrom DECIPHER OrientNucleotides
#' @importFrom DECIPHER LearnTaxa
#'
train_idtaxa <- function(x, max_group_size=10, max_iterations = 3,  allow_group_removal = TRUE, orient=FALSE,
                          get_lineage=FALSE, db = NULL, quiet = FALSE) {
  time <- Sys.time() # get time

  #Reformat to complete taxonomic hierarchy
  if(get_lineage & !is.null(db)){
    seqs <-  reformat_hierarchy(x, db=db, quiet=FALSE)
  } else if(get_lineage &  is.null(db)){
    stop("If get_lineage is TRUE, a db needs to be provided")
  } else  (seqs <- x)

  #Remove NA's
  if(length(names(seqs)[stringr::str_detect(names(seqs), ";NA;")]) > 1){
    remove <- names(seqs)[stringr::str_detect(names(seqs), ";NA;")]
    subset <- seqs[!names(seqs) %in% remove]
    message(paste0(length(seqs) - length(subset)," Sequences with NA's in taxonomy removed"))
    seqs <- subset
  }

  # Convert to DNAstringset for DECIPHER
  if(methods::is(seqs, "DNAbin")){
    seqs <-  DNAbin2DNAstringset(seqs)
  }

  # Remove gaps
  if(!quiet){ message("Removing gaps")}
  seqs <- DECIPHER::RemoveGaps(seqs)

  if(orient){
    if(!quiet){ message("Orienting sequences")}
    seqs <- DECIPHER::OrientNucleotides(seqs)
  }

  # As taxonomies are encoded in the sequence names rather than a separate file, use:
  taxid <- NULL

  #Add Root Rank
  names(seqs) <- names(seqs)  %>% stringr::str_replace(";[;]*", ";Root;")

  # obtain the taxonomic assignments
  groups <- names(seqs) # sequence names

  # assume the taxonomy begins with 'Root;'
  groups <- gsub("(.*)(Root;)", "\\2", groups) # extract the group label
  groupCounts <- table(groups)
  u_groups <- names(groupCounts) # unique groups

  if (!quiet) {message(paste0(length(u_groups), " unique species"))}

  # Pruning training set
  remove <- logical(length(seqs))
  for (i in which(groupCounts > max_group_size)) {
    index <- which(groups==u_groups[i])
    keep <- sample(length(index),
                   max_group_size)
    remove[index[-keep]] <- TRUE
  }
  if (!quiet) {message(paste0(sum(remove), " sequences pruned from over-represented groups"))}

  # Training the classifier
  probSeqsPrev <- integer() # suspected problem sequences from prior iteration

  train <- seqs[!remove]
  taxonomy <- gsub("(.*)(Root;)", "\\2", names(train))
  remove <- logical(length(train))

  for (i in seq_len(max_iterations)) {
    if (!quiet) {message("Training iteration: ", i, "\n", sep="")}
    # train the classifier
    trainingSet <- DECIPHER::LearnTaxa(train[!remove], taxonomy[!remove], taxid)

    # look for problem sequences
    probSeqs <- trainingSet$problemSequences$Index
    if (length(probSeqs)==0) {
      if (!quiet) {message("No problem sequences remaining.\n")}
      break
    } else if (length(probSeqs)==length(probSeqsPrev) &&
               all(probSeqsPrev==probSeqs)) {
      if (!quiet) {message("Iterations converged.\n")}
      break
    }
    if (i==max_iterations)
      break
    probSeqsPrev <- probSeqs

    # remove any problem sequences
    index <- which(!remove)[probSeqs]
    remove[index] <- TRUE # remove all problem sequences
    if (!allow_group_removal) {
      # replace any removed groups
      missing <- !(u_groups %in% groups[!remove])
      missing <- u_groups[missing]
      if (length(missing) > 0) {
        index <- index[groups[index] %in% missing]
        remove[index] <- FALSE # don't remove
      }
    }
  }
  if (!quiet) {message(paste0(sum(remove), " sequences removed"))}
  if (!quiet) {message(paste0(length(probSeqs), " problem sequences remaining"))}

  time <- Sys.time() - time
  if (!quiet) (message(paste0("trained IDTAXA on ", length(x), " sequences in ", format(time, digits = 2))))
  return(trainingSet)
}


# taxonomy_to_newick ------------------------------------------------------
#' tax2tree
#'
#' @param x a DNAbin object or an object coercible to DNAbin
#' @param ranks The taxonomic ranks currently assigned to the names
#' @param depth The depth within the tree to return. Default NULL to return lowest input rank
#' @param summarise Select a taxonomic rank to summarise below. Default is FALSE to return a tree of the same size as input
#' @param output The output to return, options are:
#' "phylo" to return an ape phylo object. If summarise is set to a rank, the number of distinct taxa at that rank numbers are appended to the tip labels if append_sum is TRUE, or returned as an attribute and can be accessed with attr(tree, "sum")
#' "data.tree" to return a data.tree node object
#' "treedf" to return a simplified tree with taxon summaries in a data frame
#' "newick" to return a newick file for visualisation in other software. If summarise is set to a rank, the number of distinct taxa at that rank numbers are returned as an attribute and can be accessed with attr(tree, "sum")
#' @param replace_bads An option to automatically replace invalid newick characters with '-'
#' @param append_sum An option to append the summary number to the output trees tip labels when summarise is a rank and output is 'phylo'.
#' If FALSE, summary numbers are returned as an attribute and can be accessed with attr(tree, "sum")
#'
#' @return a phylo, data.tree, newick, or treedf object
#' @import stringr
#' @import dplyr
#' @importFrom tidyr unite
#' @importFrom phytools read.newick
#' @importFrom magrittr set_colnames
#' @importFrom data.tree as.Node
#' @importFrom data.tree ToNewick
#' @importFrom data.tree ToDataFrameTree
#' @importFrom ape as.DNAbin
#' @importFrom ape base.freq
#' @importFrom methods is
#' @export
#'
#'
tax2tree <- function(x, ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
                     depth=NULL, summarise = FALSE,  output="phylo", replace_bads=FALSE, append_sum = TRUE){

  # start timer
  time <- Sys.time()

  #Checks
  if (!methods::is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
    if (all(is.na(ape::base.freq(x)))) {stop("Error: Object is not coercible to DNAbin \n")}
  }
  if (!output %in% c("phylo", "treedf", "data.tree", "newick")) { stop("Output must be 'phylo', 'data.tree', 'newick' or 'treedf'")}

  if(any(stringr::str_detect(names(x) %>% stringr::str_remove("\\|.*$"), "\\,|\\;|\\:|\\(|\\)|\\[|\\]")) && output %in% c("newick", "phylo")){
    if(replace_bads){
      names(x) <- names(x) %>%
        stringr::str_split_fixed("\\|", n=2) %>%
        as.data.frame() %>%
        dplyr::mutate(
          V1 = V1 %>%
            stringr::str_remove("\\|.*$") %>%
            stringr::str_replace_all("\\,|\\;|\\:|\\(|\\)|\\[|\\]", "_")
        ) %>%
        tidyr::unite("names", everything(), sep="|") %>%
        dplyr::pull(names)

    }else{
      stop("Sequence accessions contain one or more of the invalid characters , ; : ( ) [ ] charcters, set replace_bads to TRUE to replace with '-' ")
    }
  }

  # leave an autodetect for ranks?
  ranks <- stringr::str_to_lower(ranks)

  if (!is.null(depth)){
    ranks <- ranks[1:depth]
  }

  if(methods::is(summarise, "character")){
    summarise <- stringr::str_to_lower(summarise)

    if(summarise %in% ranks){
      groupranks <- ranks[1:match(summarise, ranks)]
    } else {
      stop("Summarise must be one of the values in ranks")
    }
  }
  #Get taxonomic lineage
  lineage <- names(x) %>%
    stringr::str_remove(pattern=";$") %>%
    stringr::str_split_fixed(";", n=Inf) %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("Acc", ranks ))

  lineage <- lineage[,!is.na(colnames(lineage))]

  if(methods::is(summarise, "character")){
    lineage <- lineage %>%
      dplyr::select(-Acc) %>%
      dplyr::group_by_at(groupranks) %>%
      dplyr::summarise(sum = dplyr::n())

    tipsums <- lineage %>%
      dplyr::pull(sum)
    names(tipsums) <- lineage%>%
      dplyr::pull(!!summarise)

    lineage <- lineage %>%
      tidyr::unite(col=pathString, !!groupranks, sep="/") %>%
      dplyr::mutate(pathString = paste0("Root/", pathString))%>%
      data.tree::as.Node(.)

  } else if(isFALSE(summarise)){
    lineage <- lineage %>%
      mutate(Acc = Acc %>% str_remove("\\|.*$")) %>%
      tidyr::unite(col=pathString, !!ranks, Acc, sep="/") %>%
      dplyr::mutate(pathString = paste0("Root/", pathString)) %>%
      data.tree::as.Node(.)
  } else(
    stop("Summarise must refer to a taxonomic rank, or be FALSE")
  )

  if (output=="phylo"){
    nwk <- lineage  %>%
      data.tree::ToNewick(heightAttribute = NULL)
    out <-  phytools::read.newick(textConnection(nwk))
    if(methods::is(summarise, "character") & append_sum){
      out$tip.label <- paste0(out$tip.label, "-", tipsums)
      attr(out, "sum") <- tipsums
    }else if(methods::is(summarise, "character") & !append_sum){
      attr(out, "sum") <- tipsums
    }
  } else if (output=="treedf" & methods::is(summarise, "character")){
    out <- data.tree::ToDataFrameTree(lineage, "sum")
  } else if (output=="treedf" & !methods::is(summarise, "character")){
    out <- data.tree::ToDataFrameTree(lineage)
  } else if (output=="data.tree"){
    out <- lineage
  } else if (output =="newick"){
    out <- lineage  %>%
      data.tree::ToNewick(heightAttribute = NULL)
    if(methods::is(summarise, "character")){
      attr(out, "sum") <- tipsums
    }
  }

  time <- Sys.time() - time
  message(paste0("Generated a taxonomic tree for ", length(x), " Sequences in ", format(time, digits = 2)))

  return(out)
}


# tax2phylo ---------------------------------------------------------------
#' tax2phylo
#'
#' @description A newer, faster tax2tree that grafts tips. Only returns a phylo objec
#' @param x a DNAbin object or an object coercible to DNAbin
#' @param ranks The taxonomic ranks currently assigned to the names
#' @param depth The depth within the ranks to return. Default all to return lowest input rank
#' @param summarise Summarise the number of sequences at depth.
#' @param append_sum Append the summary number to the output trees tip labels when summarise is TRUE
#' @param hex Convert the accessions to hexadecimal format to avoid newick incompatible characters. Can be reversed with hex2acc
#' @param resolve_poly Randomly resolve polytomies in tree using `ape::multi2di`. use 'all' to resolve all polytomies, or 'upper' to only resolve upper branches, useful for making constraints trees
#'
#' @return
#' @export
#' @import stringr
#' @import dplyr
#' @importFrom tidyr unite
#' @importFrom magrittr set_colnames
#' @importFrom phytools read.newick
#' @importFrom data.tree as.Node
#' @importFrom data.tree ToNewick
#' @importFrom data.tree ToDataFrameTree
#' @importFrom ape as.DNAbin
#' @importFrom ape base.freq
#' @importFrom ape multi2di
#' @importFrom methods is
#'
tax2phylo <- function(x, ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), depth="all",
                      summarise = FALSE, append_sum = TRUE, hex=FALSE, resolve_poly=FALSE){

  # start timer
  time <- Sys.time()

  #Checks
  if (!methods::is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
    if (all(is.na(ape::base.freq(x)))) {stop("Error: Object is not coercible to DNAbin \n")}
  }

  if(isTRUE(hex)){
    x <- acc2hex(x)
  }

  if(any(stringr::str_detect(names(x) %>% stringr::str_remove("\\|.*$"), "\\,|\\;|\\:|\\(|\\)|\\[|\\]"))) {
    stop("Sequence accessions contain one or more of the invalid characters , ; : ( ) [ ] charcters, set hex to TRUE to convert to hexadecimal encoding ")
  }

  if(is.character(resolve_poly) & !resolve_poly %in% c("upper", "all")){
    stop("resolve_poly must either be FALSE, 'upper', or 'all ")
  }

  # leave an autodetect for ranks?
  ranks <- stringr::str_to_lower(ranks)
  depth <- stringr::str_to_lower(depth)

  if (depth %in% ranks){
    groupranks <- ranks[1:match(depth, ranks)]
    message("Ranks included at depth ", depth," : ", paste(ranks, collapse=" "))
  } else if(depth==all){
    groupranks <- ranks
  } else {
    stop("Depth must be 'all' or one of the values in ranks")
  }

  #Get taxonomic lineage
  lineage <- names(x) %>%
    stringr::str_remove(pattern=";$") %>%
    stringr::str_split_fixed(";", n=Inf) %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("Acc", ranks ))

  #lineage <- lineage[,!is.na(colnames(lineage))]

  if(isTRUE(summarise)){
    lineage <- lineage %>%
      dplyr::select(-Acc) %>%
      dplyr::group_by_at(groupranks) %>%
      dplyr::summarise(sum = dplyr::n())

    tipsums <- lineage %>%
      dplyr::pull(sum)
    names(tipsums) <- lineage%>%
      dplyr::pull(!!depth)

    tree <- lineage %>%
      tidyr::unite(col=pathString, !!groupranks, sep="/") %>%
      dplyr::mutate(pathString = paste0("Root/", pathString))%>%
      data.tree::as.Node(.) %>%
      data.tree::ToNewick(heightAttribute = NULL) %>%
      textConnection() %>%
      phytools::read.newick()

    #Append summaries of lower ranks
    if(append_sum){
      tree$tip.label <- paste0(tree$tip.label, "-", tipsums)
      attr(tree, "sum") <- tipsums
    }else if(!append_sum){
      attr(tree, "sum") <- tipsums
    }

    #Resolve polytomies in tree
    if(resolve_poly == "upper"){
      tree <- ape::multi2di(tree)
    }

  } else if(isFALSE(summarise)){
    # get mapping table for grafting
    mapping <- lineage %>%
      dplyr::mutate(Acc = Acc %>% str_remove("\\|.*$")) %>%
      dplyr::select(tail(groupranks,1), Acc)%>%
      mutate_all(~str_replace_all(.x," ", "_"))

    upper_tree <- lineage %>%
      dplyr::select(!!groupranks) %>%
      dplyr::distinct() %>%
      tidyr::unite(col=pathString, !!groupranks, sep="/") %>%
      dplyr::mutate(pathString = paste0("Root/", pathString)) %>%
      data.tree::as.Node(.,
                         table,
                         pathName = "pathString",
                         pathDelimiter = "/",
                         colLevels = NULL,
                         na.rm = TRUE,
                         check = "check")  %>%
      data.tree::ToNewick(heightAttribute = NULL) %>%
      textConnection() %>%
      phytools::read.newick()

    #Resolve polytomies in upper tree
    if(resolve_poly == "upper"){
      upper_tree <- ape::multi2di(upper_tree)
    }

    # Graft tips onto upper tree
    tree <- graft_tips(upper_tree, mapping)
  }else(
    stop("Summarise must be TRUE or FALSE")
  )
  #Resolve polytomies in tree
  if(resolve_poly == "all"){
    tree <- ape::multi2di(tree)
  }
  time <- Sys.time() - time
  message(paste0("Generated a taxonomic tree for ", length(x), " Sequences in ", format(time, digits = 2)))
  return(tree)
}


# graft_tips --------------------------------------------------------------
#' Graft new tips onto an existing tree
#'
#' @param tree A phylo object
#' @param mapping A mapping data frame containing two columns.
#' The left column must be identical to the tip labels of the phylogeny, the right column contains new tips to be grafted.
#' @param branch_length A numeric giving the branch lengths for the new polytomic tips that are added.
#'
#' @return
#' @export
#' @importFrom ape read.tree
#' @importFrom ape compute.brlen
#' @importFrom ape bind.tree
#'
#'
graft_tips <- function(tree, mapping, branch_length = 0){
  if (!inherits(tree, "phylo")) stop("'tree' is not of class 'phylo'")

  info <- apply(mapping, 2, setequal, y = tree$tip.label)
  if (!any(info)) stop("'mapping' is not congruent with tiplabels of 'tree'")

  #Split into current tips to be grafted
  mapping <- split(mapping[, !info], mapping[, info])
  add_tip <- function(z){
    z <- ape::read.tree(text = paste0("(", paste(z, collapse = ","), ");"))
    ape::compute.brlen(z, branch_length) ## set branch lengths to branch_length
  }
  combs <- lapply(mapping, add_tip)
  for (i in seq_along(combs)){
    tree <- ape::bind.tree(tree, combs[[i]], which(tree$tip.label == names(combs)[i]))
  }
  return(tree)
}


# Filter sequences by taxonomy --------------------------------------------
#' Filter sequences by taxonomy
#'
#' @param x A DNAbin or DNAString with heirarchial taxonomy
#' @param filt_rank The taxonomic rank to subset at i.e. 'Class'
#' @param filt_value A taxonomic name or vector of taxonomic names to subset filt_rank at at i.e. c('Insecta', 'Arachnida')
#' @param ranks The taxonomic ranks currently assigned to the sequences
#'
#' @return
#' @export
#' @import stringr
#' @importFrom tidyr unite
#' @import dplyr
#' @importFrom magrittr set_colnames
#'
#'
filter_by_tax <- function(x, filt_rank, filt_value, ranks=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")){

  # Convert to DNAbin
  if (!class(x) %in% c("DNAbin", "DNAStringSet")  ) {
    stop("x must be a DNAbin or DNAstringset \n")
  }

  ranklength <- length(ranks)
  if(!stringr::str_count(names(x)[[1]] %>% stringr::str_remove(";$"), ";")== ranklength){
    stop("Error, the number of semicolon delimiters do not match the length of ranks")
  }

  filtnames <- names(x) %>%
    stringr::str_remove(";$") %>%
    stringr::str_split_fixed(";", n=(ranklength+1)) %>%
    as.data.frame(stringsAsFactors=FALSE) %>%
    magrittr::set_colnames(c("acc", tolower(ranks))) %>%
    dplyr::filter(!!sym(tolower(filt_rank)) %in% filt_value)%>%
    tidyr::unite(col="names", dplyr::everything(), sep=";") %>%
    pull(names)

  out <- x[names(x) %in% filtnames]

  if(length(filt_value) > 10){
    out_message <- paste0(paste(filt_value[1:10], collapse= ", "),"...[TRUNCATED]")
  } else(
    out_message <- paste(filt_value, collapse= ", ")
  )

  message(length(x) - length(out), " sequences which did not pass the filter '", filt_rank, " %in% ", out_message, "' removed")
  return(out)
}

# LCA probs ---------------------------------------------------------------

#' Probabilitys of sharing a rank as a function of sequence identity
#'
#' @param x a DNAbin object or an object coercible to DNAbin
#' @param method The distance matrix computation method to use,
#' accepts "mbed" which computes a distance matrix from each sequence to a subset of 'seed' sequences using the method outlined in Blacksheilds et al (2010).
#' This scales well to big datasets, alternatively "kdist" computes the full n * n distance matrix.
#' @param k integer giving the k-mer size used to generate the input matrix for k-means clustering.
#' @param nstart value passed to  nstart passed to kmeans. Higher increases computation time but can improve clustering accuracy considerably.
#' @param ranks The taxonomic ranks currently assigned to the names
#' @param delim The delimiter used between ranks
#'
#' @return
#' @export
#' @import stringr
#' @import dplyr
#' @importFrom kmer kdistance
#' @importFrom kmer mbed
#' @importFrom tidyr separate
#' @importFrom ape as.DNAbin
#' @importFrom utils combn
#'
lca_probs <- function(x, method="mbed",  k=5, nstart = 20, ranks=c("kingdom", "phylum", "class", "order", "family", "genus", "species"), delim=";"){
  # Convert to DNAbin
  if (!methods::is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
  }

  # Replace any trailing delimiters
  names(x) <- names(x) %>%
    stringr::str_remove(";$")
  # Count number of remaining delimiters
  ndelim <- stringr::str_count(names(x)[1], ";")

  #check if delims match ranks
  if(any(ndelim <=1)) {
    stop("Error: Needs to be in hierarchial format first, run get_ncbi_lineage or get_ott_lineage")
  } else if (!ndelim==length(ranks)){
    stop("Error: Number of delimiters does not match number of ranks")
  }

  # Create distance matrix - could probably do this using random samples
  if(method == "kdist"){
    dist <- as.matrix(kmer::kdistance(x, k=k, method="edgar", residues="DNA"))
  } else if(method == "mbed"){
    dist <- as.matrix(kmer::mbed(x, k=k, residues="DNA")[,])
  }

  #convert to pairwise distances
  xy <- t(utils::combn(colnames(dist), 2))
  pw <- data.frame(xy, dist=(100 - round(dist[xy] * 100)))

  # subset to the values in loop

  sim <- sort(unique(pw$dist), decreasing = TRUE)
  simlist <- vector("list", length=length(sim))
  s=1
  for (s in 1:length(sim)) {

    subsets <- pw %>%
      filter(dist==sim[s])

    df1 <- subsets %>%
      tidyr::separate(X1, into=c("Acc",ranks), sep=";") %>%
      dplyr::select(rev(ranks))

    df2 <- subsets %>%
      tidyr::separate(X2, into=c("Acc",ranks), sep=";") %>%
      dplyr::select(rev(ranks))

    #Get all shared ranks
    logidf <- as.data.frame(df1 == df2)

    #Get lowest common rank
    keepvec <- apply(logidf, 1, which.max)
    rows <- seq(logidf[,1])
    #cols <- seq(logidf[1,])
    #unselect <- matrix(ncol=2,c(rep(rows, length(cols)), sort(rep(cols, length(rows)))))
    select <- matrix(ncol=2,c(rows, keepvec))

    logidf[select] <- "KEEP"
    logidf[!logidf=="KEEP"] <- 0
    logidf[logidf=="KEEP"] <- 1
    logidf <- logidf  %>%
      mutate_all(as.numeric) %>%
      colSums() / length(rows)

    simlist[[s]]  <- data.frame(rank=names(logidf), prob=logidf, sim=sim[s])

  }
  names(simlist) <- sim
  out <- dplyr::bind_rows(simlist) %>%
    dplyr::group_by(rank, sim) %>%
    dplyr::summarise(prob = mean(prob))
  return(out)
}


# Accessions from fasta ---------------------------------------------------
#' Genbank accessions from fasta
#'
#' @param x A filepath or vector of filepaths to fasta files
#'
#' @return
#' @export
#' @import dplyr
#' @importFrom Biostrings fasta.index
#'
acc_from_fasta <- function(x) {
  if(!all(file.exists(x))){
    stop("Input must be a filepath or vector of filepaths to fasta files")
  }
  out <- Biostrings::fasta.index(x) %>%
    dplyr::mutate(acc = desc %>%
                    stringr::str_remove(pattern="\\|.*$")) %>%
    dplyr::pull(acc)

  return(out)
}



# Function to generate random sequences
#' Generate random sequences
#'
#' @param n The number of sequences to generate
#' @param length The length of the generated sequences
#' @param alphabet The DNA alphabet to draw from, default is A,G,C,T
#'
#' @return
#' @export
#'
random_seq <- function(n, length, alphabet = c("A","G","T","C")){
  out <- seq(1, n, 1) %>%
    purrr::map2(length, function(x,y ){
      paste(sample(alphabet, y, replace = T), collapse="")
    }) %>%
    char2DNAbin()
  names(out) <- make.unique(rep("Seq", n), sep="_")
  return(out)
}

# Accession to hexadecimal coding ----------------------------------------------
#' Convert accession number to hexadecimal coding
#' @description Hexadecimal coding is useful for ensuring that accession numbers dont contain characters that are considered illegal for newick strings
#'
#' @param x A DNAbin or DNAStringSet with names formatted Accession|taxid;taxonomy.
#' Or a character vector with names formatted Accession|taxid;taxonomy
#' Or a character vector of accessions
#' @param force override checks if string is already hexadecimal
#'
#' @return
#' @export
#' @import stringr
#' @importFrom methods is
#'
acc2hex <- function(x, force=FALSE){
  #Check input format
  if (methods::is(x, "DNAbin") | methods::is(x, "DNAStringSet")) {
    message("Input is ", class(x),", assuming header is in 'Accession|taxid;taxonomy' format")
    #Check if already hexed
    if(.ishex(names(x)[[1]] %>% stringr::str_remove("\\|.*$")) & isFALSE(force)){
      stop("Accessions are  already in hexadecimal format, set force = TRUE to override this error")
      }
    acc <- names(x) %>%
      stringr::str_remove("\\|.*$") %>%
      purrr::map_chr(.hex)
    names(x) <- names(x) %>%
      stringr::str_remove("^.*(?=\\|)") %>%
      paste0(acc, .)
  }else  if (methods::is(x, "character") && (stringr::str_detect(x, "\\|") & stringr::str_detect(x, ";"))) {
    message("Detected | and ; delimiters, assuming 'Accession|taxid;taxonomy' format")
    if(.ishex(x[[1]] %>% stringr::str_remove("\\|.*$")) & isFALSE(force)){
      stop("Accessions are  already in hexadecimal format, set force = TRUE to override this error")
    }
    acc <- x %>%
      stringr::str_remove("\\|.*$") %>%
      purrr::map_chr(.hex)
    x <- x %>%
      stringr::str_remove("^.*(?=\\|)") %>%
      paste0(acc, .)
  } else  if (methods::is(x, "character") && !(stringr::str_detect(x, "\\|") && stringr::str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of accessions")
    if(.ishex(x[[1]]) & isFALSE(force)){
      stop("Accessions are  already in hexadecimal format, set force = TRUE to override this error")
    }
    x <- x %>%
      purrr::map_chr(.hex)
  }else (stop("x must be DNA bin, DNAStringSet or vector of accesssions"))
  return(x)
}

# Hexadecimal coding to accession -----------------------------------------
#' Hexadecimal coding to accession numberD
#' @description Hexadecimal coding is useful for ensuring that accession numbers dont contain characters that are considered illegal for newick strings
#' @param x A DNAbin or DNAStringSet with names formatted Accession|taxid;taxonomy.
#' Or a character vector with names formatted Accession|taxid;taxonomy
#' Or a character vector of accessions
#' @param force override checks if string is already hexadecimal
#' @return
#' @export
#' @import stringr
#' @importFrom methods is
#'
hex2acc <- function(x, force=FALSE){
  #Check input format
  if (methods::is(x, "DNAbin") | methods::is(x, "DNAStringSet")) {
    message("Input is ", class(x),", assuming header is in 'Accession|taxid;taxonomy' format")
    if(!.ishex(names(x)[[1]] %>% stringr::str_remove("\\|.*$")) & isFALSE(force)){
      stop("Accessions are already in alphanumeric, set force = TRUE to override this error")
    }
    acc <- names(x) %>%
      stringr::str_remove("\\|.*$") %>%
      purrr::map_chr(.unhex)
    names(x) <- names(x) %>%
      str_remove("^.*(?=\\|)") %>%
      paste0(acc, .)
  }else  if (methods::is(x, "character") && (stringr::str_detect(x, "\\|") & stringr::str_detect(x, ";"))) {
    message("Detected | and ; delimiters, assuming 'Accession|taxid;taxonomy' format")
    if(!.ishex(x[[1]] %>% stringr::str_remove("\\|.*$")) & isFALSE(force)){
      stop("Accessions are already in alphanumeric, set force = TRUE to override this error")
    }
    acc <- x %>%
      stringr::str_remove("\\|.*$") %>%
      purrr::map_chr(.unhex)
    x <- x %>%
      stringr::str_remove("^.*(?=\\|)") %>%
      paste0(acc, .)
  } else  if (methods::is(x, "character") && !(stringr::str_detect(x, "\\|") && stringr::str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of accessions")
    if(!.ishex(x[[1]]) & isFALSE(force)){
      stop("Accessions are already in alphanumeric, set force = TRUE to override this error")
    }
    x <- x %>%
      purrr::map_chr(.unhex)
  }else (stop("x must be DNA bin, DNAStringSet or vector of accesssions"))
  return(x)
}


# hexadecimal translation -------------------------------------------------
#' convert alphanumerics to hexadecimal
#' internal taxreturn function
#' @param y a character string
#'
#' @return
#'
#'
.hex <- function(y){
  paste(charToRaw(y),  collapse="")
}

#' convert hexadecimal strings to alphanumeric
#' internal taxreturn function
#' @param y a character string
#'
#' @return
#'
#'
.unhex <- function(y){
  h <- sapply(seq(1, nchar(y), by=2), function(x) substr(y, x, x+1))
  rawToChar(as.raw(strtoi(h, 16L)))
}


#' Check whether a string is hexadecimal
#' internal taxreturn function
#' @param y a character string
#'
#' @return
#'
#'
.ishex <- function(y) {
  h <- sapply(seq(1, nchar(y), by = 2), function(x) substr(y, x, x + 1))
  if(any(is.na(strtoi(h, 16L)))){
    return(FALSE)
  } else{
    return(TRUE)
  }
}


# Multithread -------------------------------------------------------------
#' Setup multithreading
#'
#' @param multithread Number of cores
#' @param quiet Whether progress should be printed to console
#'
#' @return
#' @import future
#'
#' @examples
setup_multithread <- function(multithread, quiet=FALSE){
  ncores <- future::availableCores()
  if(isTRUE(multithread)){
    cores <- ncores-1
    if(!quiet){message("Multithreading with ", cores, " cores")}
    future::plan(future::multisession, workers=cores)
  } else if (is.numeric(multithread) & multithread > 1){
    cores <- multithread
    if(cores > ncores){
      cores <- ncores
      warning("The value provided to multithread is higher than the number of cores, using ", cores, " cores instead")
    }
    if(!quiet){message("Multithreading with ", cores, " cores")}
    future::plan(future::multisession, workers=cores)
  } else if(isFALSE(multithread) | multithread==1){
    future::plan(future::sequential)
  } else (
    stop("Multithread must be a logical or numeric vector of the numbers of cores to use")
  )
}

#' Setup parallel
#'
#' @param multithread Number of cores
#' @param quiet Whether progress should be printed to console
#'
#' @return
#' @importFrom future availableCores
#'
#' @examples
setup_para <- function(multithread, quiet=FALSE){
  ncores <- future::availableCores() -1
  if(isTRUE(multithread)){
    cores <- ncores
    if(!quiet){message("Multithreading with ", cores, " cores")}
  } else if (is.numeric(multithread) & multithread > 1){
    cores <- multithread
    if(cores > ncores){
      cores <- ncores
      warning("The value provided to multithread is higher than the number of cores, using ", cores, " cores instead")
    }
    if(!quiet){message("Multithreading with ", cores, " cores")}
  } else if(isFALSE(multithread) | multithread==1){
    cores <- 1
  } else (
    stop("Multithread must be a logical or numeric vector of the numbers of cores to use")
  )
  return(cores)
}

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL