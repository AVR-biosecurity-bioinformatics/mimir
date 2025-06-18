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

# replace each space in sequence headers of fasta with the string "!?!?"
sed '/^>/ s/ /!?!?/g' $4 > renamed.fasta

# create nucleotide DB
mmseqs createdb \
	renamed.fasta \
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
	clusters_pre.tsv \
	--threads $2 

# convert tsv back to original sequence names
sed 's/!?!?/ /g' clusters_pre.tsv > clusters.tsv

rm clusters_pre.tsv
rm renamed.fasta