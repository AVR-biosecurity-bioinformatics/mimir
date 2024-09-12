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
