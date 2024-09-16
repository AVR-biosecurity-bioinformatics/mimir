# dev notes


run pipeline
```
module load Java/17
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds
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
- add dynamic resources 


Less important
- add 'chunk_size' parameter to chunk .fasta files into this many sequences after fetching