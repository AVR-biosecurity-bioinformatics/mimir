# dev notes


run pipeline
```
module load Java/17
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds

# entrez key (example only)
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxon 7215 --chunk_rank species

# testing whole-order (Archaeognatha)
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxon 29994 --chunk_rank family


```


### TODO

Important
- make process that creates a PHMM from a given .fasta of marker sequences
- implement deduplication process to check for same accessions in BOLD and NCBI 
- add pipeline parameters that control process function parameters (ie. map_to_model)
- add pipeline parameter for marker search query 
- add pipeline parameter for taxon search parameter 
- add ability to add internal sequences through .fasta file
    - what format does the name of each sequence need to be in?
    - preferentially retain sequences in PRUNE_GROUPS step
- add dynamic resources based on file size/number of 
- check that chunked alignments (to PHMM) are equivalent to alignments on combined .fasta
- query ncbi taxdump db for all ranks at certain level, then all ranks not captured etc. -- to produce taxa chunks for fetching sequences


Less important
- add 'chunk_size' parameter to chunk .fasta files into this many sequences after fetching
- add process that tracks the number of sequences at each step of the pipeline
- add option to retrieve a small number of outgroup taxa