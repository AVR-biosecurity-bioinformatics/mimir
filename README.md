# `mimir`: a DNA barcode reference database curation pipeline

`mimir` is a [Nextflow](https://www.nextflow.io/docs/latest/index.html)-based pipeline for reference database curation of DNA barcode sequences, such as those required for metabarcoding experiments. It is partially a successor to [`taxreturn`](https://github.com/alexpiper/taxreturn).

This pipeline is being developed by a team at [Agriculture Victoria Research](https://agriculture.vic.gov.au/), as a part of the [National Grains Diagnostic & Surveillance Initiative (NGDSI)](https://grdc.com.au/grdc-investments/investments/investment?code=DEE2305-004RTX). 

`mimir` is a sibling pipeline to [`freyr`](https://github.com/AVR-biosecurity-bioinformatics/freyr), which handles sequencing data analysis for metabarcoding experiments. 

> This pipeline is currently **UNFINISHED** and being actively developed, with no guarantee that the code is stable or usable!


### Typical commands

Fetching COI sequences belonging to the genus *Drosophila* (NCBI taxid: 7215) from GenBank only:

```
NXF_VER=23.05.0-edge nextflow run . \
    --target_taxon 7215 \
    --target_rank genus \
    --marker COI \
    --use_bold false
```

Fetching COI sequences from the same taxon using GenBank and BOLD (using a manually generated URL to download BOLD database):

```
NXF_VER=23.05.0-edge nextflow run . \
    --target_taxon 7215 \
    --target_rank genus \
    --marker COI \
    --bold_db_url "https://www.boldsystems.org/index.php/API_Datapackage?id=BOLD_Public.06-Sep-2024&uid=166f4b93030986"
```

Fetching COI sequences from GenBank and combining with curated 'internal' sequences in the file `my_sequences.fasta`:

```
NXF_VER=23.05.0-edge nextflow run . \
    --target_taxon 7215 \
    --target_rank genus \
    --marker COI \
    --use_bold false \
    --internal_seqs my_sequences.fasta
```

Fetching and trimming COI sequences to a primer pair, keeping sequences if >=200 bp after trimming:

```
NXF_VER=23.05.0-edge nextflow run . \
    --target_taxon 7215 \
    --target_rank genus \
    --marker COI \
    --use_bold false \
    --trim_to_primers \
    --primer_fwd GGDACWGGWTGAACWGTWTAYCCHCC \
    --primer_rev GTRATWGCHCCDGCTARWACWGG \
    --min_length_trimmed 200
```

Fetching COI sequences from a list of taxa in `target_list_file.csv`:

```
NXF_VER=23.05.0-edge nextflow run . \
    --target_list target_list_file.csv \
    --marker COI \
    --use_bold false
```

Keeping only sequences that are classified at every key rank:
```
NXF_VER=23.05.0-edge nextflow run . \
    --target_taxon 7215 \
    --target_rank genus \
    --marker COI \
    --use_bold false \
    --remove_unclassified any_ranks
```





### Notes on usage

#### Fetching sequences from GenBank

It is strongly recommended users obtain an [NCBI API key](https://support.nlm.nih.gov/kbArticle/?pn=KA-05317) and supply it to the pipeline using `--entrez_key`. This will increase the rate at which sequences can be fetched from GenBank and will speed up the pipeline. 


#### Fetching sequences from BOLD 

Best way to access Barcode of Life Data ([BOLD](https://www.boldsystems.org/index.php)) data at the moment is with a direct link to the latest [complete database package](https://bench.boldsystems.org/index.php/datapackages/Latest). You'll need to [create an account](https://bench.boldsystems.org/index.php/MAS_Management_NewUserApp) with BOLD, which will allow you to generate a URL that will download the complete database package (>2 GB compressed) for 24 hours. Input an active URL into the pipeline in the following format, making sure to wrap the URL in quotes: `--bold_db_url "https://www.boldsystems.org/index.php/API_Datapackage?id=BOLD_Public.06-Sep-2024&uid=166f4b93030986"`

If you have already downloaded a complete BOLD data package, either from a previous pipeline run or a manual download, you can supply a path to either the compressed `.tar.gz` file itself or a directory where the compressed or uncompressed file(s) are found, using the `--bold_db_path` flag. 

> TODO: Add ability for BOLD login details to be given to pipeline instead of generated URL.  

We are also looking into ways of downloading sequence data from BOLD using the [`bold`](https://docs.ropensci.org/bold/index.html) R package, but this method is inherently less reliable due to unpredictable rate-limiting from the BOLD website. 

#### Internal sequences

A common goal of building a reference database is to combine sequences from public databases ('external' sequences) with high-quality, curated sequences the user generated themselves or otherwise trusts ('internal' sequences). You can supply internal sequences (in a `.fasta` file) to Mimir with `--internal_seqs`, and they will be processed alongside external sequences. 

The only limitation currently is that the `.fasta` headers of the internal sequences must be in the following format: `acc|taxid;kingdom;phylum;class;order;family;genus;species`, where:
- `acc` is a unique sequence accession/id, eg. `GBMH0095-06`
- `taxid` is a numerical taxonomic ID, either from NCBI (and prefixed with `NCBI:`), BOLD (prefix: `BOLD:`) or internal use (prefix: `INTERNAL:`). For example: a sequence from *Drosophila melanogaster* would have a taxid of `NCBI:7227`
- `kingdom` to `species` is the taxonomic lineage, using names from NCBI Taxonomy. For example, a *Drosophila melanogaster* sequence would have the following lineage: `Metazoa;Arthropoda;Insecta;Diptera;Drosophilidae;Drosophila;Drosophila melanogaster`



#### Target taxa

There are two ways to tell Mimir which taxonomic group(s) (ie. targets) you would like the database to focus on. For a single taxon, use `--target_taxon` to specify a valid NCBI taxon name or UID, and `--target_rank` to specify its taxonomic rank. For multiple taxa, use `--target_list` to supply a two-column .csv that has the NCBI taxon name or UID (column 1) and the taxon rank (column 2), with one row for each taxon.

For the best compatibility between source databases (ie. Genbank and BOLD), only taxa from the following ranks should be used as targets: kingdom, phylum, class, order, family, genus, or species. This is because intermediate ranks, like "subfamily" and "tribe", are not available for all sequences. If a database comprising only an intermediate rank is needed, we recommend filtering a larger database using a known list of component taxa, or using such a list as the input to the pipeline. 

Currently, the pipeline parses input taxa through [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy), so each target taxon name must be a recognised NCBI taxon name (eg. `Insecta`), or even better, an NCBI taxid number/UID (eg. `50557`). Taxids are preferred as they are unambiguous, while taxon names are sometimes shared between different taxa (a good example is `Drosophila` belonging to three taxids: `7215`, `32281` and `2081351`). To assist with (but not completely solve) taxon name disambiguation, `--target_rank` must also be used to supply the taxonomic rank of the taxon name (eg. `class` for `Insecta`); this also helps taxonomic harmonisation between the NCBI and BOLD taxonomies. 

#### Markers/barcodes

Specify your barcode marker using `--marker`, eg. `--marker COI`. At the moment, only one marker is supported: `COI` (cytochrome oxidase 1). More markers will be available soon. The pipeline uses protein-based profile hidden Markov models (PHMMs) from [InterPro](https://www.ebi.ac.uk/interpro/) to detect, filter and trim input sequences to the specific barcode. Support for user-defined PHMMs will come in the future. 

#### Primer-based trimming

Often users will want to trim database sequences to a particular region amplified by a primer pair. Mimir does this by using user-specified primer sequences to trim the barcode's PHMM and using that to trim each sequence. To enable primer-based trimming, use `--trim_to_primers` and specify your primer sequences using `--primer_fwd` and `--primer_rev`. Base ambiguity/degeneracy is supported, and the pipeline will try to determine the correct strand orientation of the primers relative to the barcode PHMM. 