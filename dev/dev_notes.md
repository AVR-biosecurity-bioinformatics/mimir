# dev notes


run pipeline
```
NXF_VER=23.05-edge

module load Java/17
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds

# entrez key (example only; Apis + subgenus)
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 83321 --chunk_rank genus

# testing whole-order (Archaeognatha)
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 29994 --chunk_rank family

# quick test across a small order (Neuroptera)
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7516 --chunk_rank family

# testing handling of chunk_rank above target taxon rank (Drosophila melanogaster vs family)
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7227 --chunk_rank family

# Drosophilinae chunked to genus
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 43845 --chunk_rank genus

# all insects as a test
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 50557 --chunk_rank family

# testing BOLD DB downloading
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 50557 --chunk_rank family --bold_db_url "https://www.boldsystems.org/index.php/API_Datapackage?id=BOLD_Public.06-Sep-2024&uid=166f9f24c69328"

# testing BOLD DB pre-uploaded (works with both path and files)
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 50557 --chunk_rank family --bold_db_path ./input

# testing new PARSE_TARGETS on valid BOLD taxon
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Diptera --target_ranks order --bold_db_path ./input

# testing PARSE_TARGETS on invalid BOLD taxon
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Bombycoidea --target_ranks superfamily --bold_db_path ./input

# testing BOLD code
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Bombycoidea --target_ranks superfamily --bold_db_path ./input

# quick test
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Neuroptera --target_ranks order --bold_db_path ./input

# quick test 2 (Drosophila)
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7215 --target_ranks genus --bold_db_path ./input

# quick test 3 (Apis)
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7459 --target_ranks genus --bold_db_path ./input

# quick test 4 (Embioptera)
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 50657 --target_ranks order --bold_db_path ./input

# bactrocera test of clustering
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 27456 --target_ranks genus --bold_db_path ./input --cluster_rank species --cluster_threshold 0.99

# Drosophilidae test of clustering
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7214 --target_ranks family --bold_db_path ./input --cluster_rank genus --cluster_threshold 0.97 --marker COI

# Orthoptera test
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 6993 --target_ranks order --bold_db_path ./input --cluster_rank genus --cluster_threshold 0.97

# Acrididae test for NA taxids in input sequences (example KU184830.1)
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7002 --target_ranks family --bold_db_path ./input --cluster_rank genus --cluster_threshold 0.97

### testing markers
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Neuroptera --target_ranks order --bold_db_path ./input --marker COI


### testing internal sequences (correct formatting)

# Drosophila melanogaster + internal
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7227 --target_ranks species --bold_db_path ./input --marker COI --internal_seqs assets/internal_fake.fasta --cluster_rank species --cluster_threshold 0.99

# test alignment on small 
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7227 --target_ranks species --bold_db_path ./input --marker COI --internal_seqs assets/internal_fake.fasta --add_root --aligned_output

# test alignment on Neuroptera 
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Neuroptera --target_ranks order --bold_db_path ./input --marker COI --internal_seqs assets/internal_fake.fasta --add_root --aligned_output

## small tests

# Neuroptera (order), bold only
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Neuroptera --target_ranks order --bold_db_path ./input --marker COI --add_root --use_genbank false --train_idtaxa

# Drosophila melanogaster + internal, bold only
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7227 --target_ranks species --bold_db_path ./input --marker COI --internal_seqs assets/internal_fake.fasta --use_genbank false

# Aphis (genus), bold only
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 80764 --target_ranks genus --bold_db_path ./input --marker COI --use_genbank false

# Schizaphis (genus), bold only
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 13261 --target_ranks genus --bold_db_path ./input --marker COI --use_genbank false 

## big tests

# all aphids for aphid prey data (retain unclassified)
NXF_VER=23.05.0-edge nextflow run . \
	-profile basc_slurm,debug \
	--phmm_model assets/folmer_fullength_model.rds \
	--entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 \
	--target_taxa 33385 \
	--target_ranks superfamily \
    --bold_db_path ./input \
    --min_length 100 \
	--max_length 1500 \
	--marker COI \
	--add_root \
	--train_idtaxa

# all aphids for aphid prey data (remove unclassified)
NXF_VER=23.05.0-edge nextflow run . \
	-profile basc_slurm,debug \
	--phmm_model assets/folmer_fullength_model.rds \
	--entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 \
	--target_taxa 33385 \
	--target_ranks superfamily \
    --bold_db_path ./input \
    --min_length 100 \
	--max_length 1500 \
	--marker COI \
	--add_root \
	--train_idtaxa \
    --remove_unclassified any_ranks

# lepidoptera test, no internal (1 day requested) - 2111104 input
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Lepidoptera --target_ranks order --bold_db_path ./input --marker COI --add_root 

# hemiptera test, no internal (1 day requested) - 741137 input - used 23.6G memory for REMOV_CONTAM for 696767 seqs
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Hemiptera --target_ranks order --bold_db_path ./input --marker COI --add_root 

# hemiptera BOLD only - 554389 input
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Hemiptera --target_ranks order --bold_db_path ./input --marker COI --add_root --use_genbank false

# lepidoptera test, bold only - 1506300 input
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Lepidoptera --target_ranks order --bold_db_path ./input --marker COI --add_root --use_genbank false 

# bold only insecta - 10629784 input
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Insecta --target_ranks class --bold_db_path ./input --marker COI --add_root --use_genbank false

## memory tests for REMOVE_CONTAM

# hemiptera test, no internal (1 day requested) - 741137 input - used 23.6G memory in REMOV_CONTAM for 696767 seqs
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Hemiptera --target_ranks order --bold_db_path ./input --marker COI --add_root 

# drosophila melanogaster - used 0.440G memory in REMOV_CONTAM for 137 seqs
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7227 --target_ranks species --bold_db_path ./input --marker COI --internal_seqs assets/internal_fake.fasta --use_genbank false

# drosophila genus - used 1.105G memory in REMOV_CONTAM for 19539 seqs
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7215 --target_ranks genus --bold_db_path ./input --marker COI --internal_seqs assets/internal_fake.fasta --use_genbank false

# Drosophilidae family - used 2.89G memory in REMOV_CONTAM for 68469 seqs
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7214  --target_ranks family --bold_db_path ./input --marker COI --internal_seqs assets/internal_fake.fasta --use_genbank false

# clustalo header truncation test (genus that causes issue) 
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 299289 --target_ranks genus --bold_db_path ./input --marker COI --add_root --use_genbank false

## time tests for complete pipeline (20/12/2024)
# siphonaptera BOLD only - input 600 - 12m 
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7509 --target_ranks order --bold_db_path ./input --marker COI --add_root --use_genbank false

# hemiptera BOLD only - 554389 input - 1h44m
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Hemiptera --target_ranks order --bold_db_path ./input --marker COI --add_root --use_genbank false

# Diptera BOLD only 
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Diptera --target_ranks order --bold_db_path ./input --marker COI --add_root --use_genbank false

# insecta BOLD only, no unclassified - 10629757 input (2981685 full classification) - 1 day requested
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Insecta --target_ranks class --bold_db_path ./input --marker COI --add_root --use_genbank false --remove_unclassified any_ranks

# Ascomycota (phylum) -- fungal test!
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Ascomycota --target_ranks phylum --bold_db_path ./input --marker COI --add_root --use_genbank false --remove_unclassified any_ranks


# memory usage tests (GB) - min: 0.44, mean: 0.64, median: 0.45, max: 1.75 (completed) for Insecta pre-GC additions
#  usage for below test with Drosophilidae - min: 0.44, mean: 0.63, median: 0.59, max: 1.02 (completed) post-GC
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7214  --target_ranks family --bold_db_path ./input --marker COI --internal_seqs assets/internal_fake.fasta --use_genbank false

## trimming primers testing
# Drosophila genus
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7215 --target_ranks genus --bold_db_path ./input --marker COI --internal_seqs assets/internal_fake.fasta --use_genbank false --trim_to_primers --primer_fwd GGDACWGGWTGAACWGTWTAYCCHCC --primer_rev GTRATWGCHCCDGCTARWACWGG

# siphonaptera
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa 7509 --target_ranks order --bold_db_path ./input --marker COI --add_root --use_genbank false --trim_to_primers --primer_fwd GGDACWGGWTGAACWGTWTAYCCHCC --primer_rev GTRATWGCHCCDGCTARWACWGG

# hemiptera - 6364 from prune_groups, alignment took 30s, align_primers took 7min
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Hemiptera --target_ranks order --bold_db_path ./input --marker COI --add_root --use_genbank false --remove_unclassified any_ranks --trim_to_primers --primer_fwd GGDACWGGWTGAACWGTWTAYCCHCC --primer_rev GTRATWGCHCCDGCTARWACWGG --min_length_trimmed 200

# lepidoptera - 89203 from prune_groups, alignment took 1hr 23min, align_primers his OOM error -- trying with 8 threads took... 
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Lepidoptera --target_ranks order --bold_db_path ./input --marker COI --add_root --use_genbank false --remove_unclassified any_ranks --trim_to_primers --primer_fwd GGDACWGGWTGAACWGTWTAYCCHCC --primer_rev GTRATWGCHCCDGCTARWACWGG --min_length_trimmed 200

# neuroptera
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Neuroptera --target_ranks order --bold_db_path ./input --marker COI --add_root --use_genbank false --remove_unclassified any_ranks --trim_to_primers --primer_fwd GGDACWGGWTGAACWGTWTAYCCHCC --primer_rev GTRATWGCHCCDGCTARWACWGG --min_length_trimmed 200

## testing AM's primers on neuroptera
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Neuroptera --target_ranks order --bold_db_path ./input --marker COI --add_root --use_genbank false --remove_unclassified any_ranks --trim_to_primers --primer_fwd CCHGAYATRGCHTTYCCHCG --primer_rev TCDGGRTGNCCRAARAAYCA --min_length_trimmed 350

# on lepidoptera (also a test of total pipeline)
NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,debug --phmm_model assets/folmer_fullength_model.rds --entrez_key 364ddb16f9f8fdf6133982af89d0bd762c09 --target_taxa Lepidoptera --target_ranks order --bold_db_path ./input --marker COI --add_root --use_genbank false --remove_unclassified any_ranks --trim_to_primers --primer_fwd CCHGAYATRGCHTTYCCHCG --primer_rev TCDGGRTGNCCRAARAAYCA --min_length_trimmed 350

```


### TODO

Questions
- What do we do with subspecies-level taxonomic resolution? Replace binomial with trinomial? Truncate to binomial?
    - Currently pipeline retains specific taxids (NCBI, or BOLD if can't be exactly matched) for each sequence, while populating the lineage string with only the specified ranks (ie. kingdom through species). This means a sequence identified to subspecies will have a different taxid to one identified to the same species, but will share the same lineage string. 
- How do we handle the "identification_method" field in the BOLD datbase?
    - Possibly problematic values are "BIN Taxonomy Match [date]" and "BOLD Sequence Classifier", which use BOLD data to classify the sequence and as such are not necessarily externally validated (can result in overconfidence in assignment, GIGO)
    - Could have a pipeline parameter that removes all sequences that have an ID method containing "BIN", "barcode", "BOLD", "DNA" or "tree" (but retains if tree + morphology was used?)
- How does PRUNE_GROUPS handle unclassified sequences? ie. does it treat all unclassified in a particular genus as a group, which then gets pruned to a max size? 
    - ANSWER: prune_groups function uses taxid + lineage string to determine taxonomic grouping, so "Unclassified" sequences will not be grouped together unless they share an identical taxid
    - It's possible we might want to keep unclassified sequences regardless of grouping, eg. two family-level ID sequences might be grouped but only share 90% seq ID, so are not real duplicates of each other

Important
- change REMOVE_CONTAM code so unclassified sequences at particular rank are removed before OTU clustering, then added at the end (this should speed up clustering)
- do a retry trick with QUERY_GENBANK in case of server error
- double-check SUMMARISE_TAXA is working correctly (example is quick test 4 above gives 1 species and 91 genera)
- make process that creates a PHMM from a given .fasta of marker sequences
- add pipeline parameter for marker search query 
    - add this to COUNT_GENBANK and FETCH_GENBANK processes
- add ability to add internal sequences through .fasta file
    - what format does the name of each sequence need to be in? ("accession|taxid")
        - "taxid" can be NCBI format (eg. "NCBI:11111") if assignment is known, or internal if a new species (eg. "INT:18"); if INT, lineage information must be given in format "kingdom;phylum..."
    - preferentially retain internal sequences in PRUNE_GROUPS step
- add outgroup handling
- check that chunked alignments (to PHMM) are equivalent to alignments on combined .fasta
- create "reject" channels that capture sequences that fail the phmm, stop, exact, contam and prune filters, as combined .fasta files
- allow users to skip PHMM, STOP, EXACT, CONTAM and PRUNE steps (mainly to allow quantitative comparisons between filter settings)
- allow users to input a "fetched, unfiltered" .fasta file (ie. didn't filter anything) to allow quantitative comparisons of the filtering steps
- add parameter to remove all unclassified sequences from final database
- add versioning process that tracks version of BOLD db, versions of each container etc. 
- allow 'prune_groups' function to use internal sequence accessions as the "preferred" sequence names (use stringr matching)


Less important
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