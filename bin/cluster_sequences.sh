#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = task memory in kilobytes
# $4 = fasta_file
# $5 = params.cluster_threshold

# parse memory limit
TASK_MEMORY_KB=$3
SPLIT_MEMORY_LIMIT=$(( TASK_MEMORY_KB * 4 / 5 )) # 80% of total memory should go to mmseqs2: https://github.com/soedinglab/mmseqs2/wiki#memory-consumption

# create nucleotide DB
mmseqs createdb \
	$4 \
	DB \
	--dbtype 2

# linclust
mmseqs linclust \
	DB \
	DB.clustered \
	tmp \
	--min-seq-id $5 \
	--threads $2 \
	--kmer-per-seq 2000 \
	--split-memory-limit ${SPLIT_MEMORY_LIMIT}K

# tsv output
mmseqs createtsv \
	DB \
	DB \
	DB.clustered \
	clusters.tsv \
	--threads $2 
