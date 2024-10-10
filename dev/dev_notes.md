# dev notes


run pipeline
```
module load Java/17
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds

# entrez key (example only; Apis + subgenus)
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 83321 --chunk_rank genus

# testing whole-order (Archaeognatha)
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 29994 --chunk_rank family

# quick test across a small order (Neuroptera)
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7516 --chunk_rank family

# testing handling of chunk_rank above target taxon rank (Drosophila melanogaster vs family)
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7227 --chunk_rank family

# Drosophilinae chunked to genus
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 43845 --chunk_rank genus

# all insects as a test
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 50557 --chunk_rank family

# testing BOLD DB downloading
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 50557 --chunk_rank family --bold_db_url "https://www.boldsystems.org/index.php/API_Datapackage?id=BOLD_Public.06-Sep-2024&uid=166f9f24c69328"

# testing BOLD DB pre-uploaded (works with both path and files)
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 50557 --chunk_rank family --bold_db_path ./input

# testing new PARSE_TARGETS on valid BOLD taxon
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Diptera --target_ranks order --bold_db_path ./input

# testing PARSE_TARGETS on invalid BOLD taxon
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Bombycoidea --target_ranks superfamily --bold_db_path ./input

# testing BOLD code
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Bombycoidea --target_ranks superfamily --bold_db_path ./input

# quick test
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Neuroptera --target_ranks order --bold_db_path ./input

# quick test 2 (Drosophila)
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7215 --target_ranks genus --bold_db_path ./input

```


### TODO

Questions
- What do we do with subspecies-level taxonomic resolution? Replace binomial with trinomial? Truncate to binomial?
- How do we handle the "identification_method" field in the BOLD datbase?
    - Possibly problematic values are "BIN Taxonomy Match [date]" and "BOLD Sequence Classifier", which use BOLD data to classify the sequence and as such are not necessarily externally validated (can result in overconfidence in assignment, GIGO)
    - Could have a pipeline parameter that removes all sequences that have an ID method containing "BIN", "barcode", "BOLD", "DNA" or "tree" (but retains if tree + morphology was used?)

Important
- make process that creates a PHMM from a given .fasta of marker sequences
- add pipeline parameter for marker search query 
    - add this to COUNT_GENBANK and FETCH_GENBANK processes
- add ability to add internal sequences through .fasta file
    - what format does the name of each sequence need to be in?
    - preferentially retain internal sequences in PRUNE_GROUPS step
- modify fetch_seqs/fetch_genbank to query taxid directly rather than taxon name to eliminate ambiguity
- add outgroup handling
- check that chunked alignments (to PHMM) are equivalent to alignments on combined .fasta
- query ncbi taxdump db for all ranks at certain level, then all ranks not captured etc. -- to produce taxa chunks for fetching sequences
    - ie. replace CHUNK_TAXON code using taxize with querying the taxdump
- find out why some sequences download without a taxid attached (and thus can't have their full taxonomic lineage resolved)
- dereplicate exact sequences before alignment to PHMM (combine .fasta chunks and then remove exact duplicates -- compare taxid as well) -- add counts to table
- modify MATCH_BOLD so it acts on the individual chunk elements output from EXTRACT_BOLD, to avoid serious memory hog issues (if the chunks are merged then acted on, memory can blow out)
    - output .fasta files (as well as the output .csvs) can then easily be concatenated in a new MERGE_BOLD step written in bash which is must less memory intensive


Less important
- add process that tracks the number of sequences at each step of the pipeline
    - count sequences from intermediate .fasta files (once moved to .fasta-only implementation)
- add option to retrieve a small number of outgroup taxa
- TRAIN_IDTAXA fails if only one taxonomic group is present (ie. one species)
    - create check for this
- REMOVE_CONTAM fails when '# testing whole-order (Archaeognatha)' command above is run, possibly due to an "NA" somewhere in a name
- add check in schema that --genetic_code is within Biostrings::GENETIC_CODE_TABLE
- make public database sourcing optional (ie. can choose Genbank or BOLD or both)
- add ability to merge existing databases together and process them 
- de-align output sequences
- output .fasta as well as .csv or .tsv for easier post-pipeline taxonomic filtering
- add parameter for fasta of excluded sequences (eg. known pseudogenes, human COI, wolbachia sequences etc.)

### NOTES

- 15m 10s to do Neuroptera database with 3356 sequences across 1791 species



### shifter tests of Docker images

```
module load shifter

shifter --image=emehinovic72/edirect:latest /bin/bash

```