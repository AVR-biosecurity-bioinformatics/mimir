# dev notes


run pipeline
```
module load Java/17
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds

# entrez key (example only; Apis + subgenus)
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxon 83321 --chunk_rank genus

# testing whole-order (Archaeognatha)
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxon 29994 --chunk_rank family

# quick test across a small order (Neuroptera)
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxon 7516 --chunk_rank family

# testing handling of chunk_rank above target taxon rank (Drosophila melanogaster vs family)
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxon 7227 --chunk_rank family

# Drosophilinae chunked to genus
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxon 43845 --chunk_rank genus


```


### TODO

Important
- make process that creates a PHMM from a given .fasta of marker sequences
- implement deduplication process to check for same accessions in BOLD and NCBI 
- add pipeline parameters that control process function parameters (ie. map_to_model)
- add pipeline parameter for marker search query 
    - add this to COUNT_GENBANK and FETCH_GENBANK processes
- add pipeline parameter for taxon search parameter 
- add ability to add internal sequences through .fasta file
    - what format does the name of each sequence need to be in?
    - preferentially retain internal sequences in PRUNE_GROUPS step
- modify fetch_seqs/fetch_genbank to query taxid directly rather than taxon name to eliminate ambiguity
- add outgroup handling
- check that chunked alignments (to PHMM) are equivalent to alignments on combined .fasta
- query ncbi taxdump db for all ranks at certain level, then all ranks not captured etc. -- to produce taxa chunks for fetching sequences
- find out why some sequences download without a taxid attached (and thus can't have their full taxonomic lineage resolved)

Less important
- add 'chunk_size' parameter to chunk .fasta files into this many sequences after fetching
- add process that tracks the number of sequences at each step of the pipeline
- add option to retrieve a small number of outgroup taxa
- TRAIN_IDTAXA fails if only one taxonomic group is present (ie. one species)
    - create check for this
- REMOVE_CONTAM fails when '# testing whole-order (Archaeognatha)' command above is run, possibly due to an "NA" somewhere in a name
- add check in schema that --genetic_code is within Biostrings::GENETIC_CODE_TABLE
- make public database sourcing optional (ie. can choose Genbank or BOLD or both)
- add ability to merge existing databases together and process them 

### NOTES

- 15m 10s to do Neuroptera database with 3356 sequences across 1791 species