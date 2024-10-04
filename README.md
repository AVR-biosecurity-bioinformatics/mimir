# `mimir`: a DNA barcode reference database curation pipeline

`mimir` is a [Nextflow](https://www.nextflow.io/docs/latest/index.html)-based pipeline for reference database curation of DNA barcode sequences, such as those required for metabarcoding experiments. It is partially a successor to [`taxreturn`](https://github.com/alexpiper/taxreturn).

This pipeline is being developed by a team at [Agriculture Victoria Research](https://agriculture.vic.gov.au/), as a part of the [National Grains Diagnostic & Surveillance Initiative (NGDSI)](https://grdc.com.au/grdc-investments/investments/investment?code=DEE2305-004RTX). 

`mimir` is a sibling pipeline to [`freyr`](https://github.com/AVR-biosecurity-bioinformatics/freyr), which handles sequencing data analysis for metabarcoding experiments. 

> This pipeline is currently **UNFINISHED** and being actively developed, with no guarantee that the code is stable or usable!




### Notes on usage

#### BOLD database 

Best way to access Barcode of Life Data ([BOLD](https://www.boldsystems.org/index.php)) data at the moment is with a direct link to the latest [complete database package](https://www.boldsystems.org/index.php/datapackages/Latest). You'll need to [create an account](https://www.boldsystems.org/index.php/MAS_Management_NewUserApp) with BOLD, which will allow you to generate a URL that will download the complete database package (>2 GB compressed) for 24 hours. Input an active URL into the pipeline in the following format, making sure to wrap the URL in quotes: `--bold_db_url "https://www.boldsystems.org/index.php/API_Datapackage?id=BOLD_Public.06-Sep-2024&uid=166f4b93030986"`

If you have already downloaded a complete BOLD data package, either from a previous pipeline run or a manual download, you can supply a path to either the compressed `.tar.gz` file itself or a directory where the compressed or uncompressed file(s) are found, using the `--bold_db_path` flag. 

> TODO: Add ability for BOLD login details to be given to pipeline instead of generated URL.  

We are also looking into ways of downloading sequence data from BOLD using the [`bold`](https://docs.ropensci.org/bold/index.html) R package, but this method is inherently less reliable due to unpredictable rate-limiting from the BOLD website. 

#### Target taxon/taxa

For the best compatibility between databases, only taxa from the following ranks should be used as `--target_taxa`: kingdom, phylum, class, order, family, genus, or species. Intermediate ranks, like "subfamily" and "tribe", are not available for all sequences; we recommend filtering the output database using a known list of component taxa, or using such a list as the input to the pipeline. 

Sometimes the taxonomic names between NCBI and BOLD differ -- by default, the pipeline parses input taxa through the NCBI taxonomy, so 