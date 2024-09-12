# dev notes


run pipeline
```
module load Java/17
NXF_VER=23.04.5 nextflow run . -profile basc_slurm,debug --phmm_model ./assets/folmer_fullength_model.rds
```


### TODO

- make process that creates a PHMM from a given .fasta of marker sequences
- implement deduplication process to check for same accessions in BOLD and NCBI 
- add pipeline parameters that control process function parameters (ie. map_to_model)