#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = db_tsv_file
# $4 = db_meta_file
# $5 = bold_targets
# $6 = marker
# $7 = min_length_input
# $8 = max_length_input

BOLD_TSV=$3
MARKER=$6
MIN_LEN=$7
MAX_LEN=$8

# remove header from bold_targets
tail -n+2 $5 > taxa.csv

## filtering steps are:
# keep rows that match marker code 
# keep rows that match the taxon name in the particular taxon rank column
# keep rows with a nucleotide sequence between min and max lengths
# keep only some columns: process id, tax id, taxonomic ranks, sequence, id method, external accession ID, and sequence run site 

# get index locations of column names for awk
head -n1 $3 | tr '\t' '\n' > fields_nl.txt # header as one-column-per-line file
MARKER_CODE_ID=$(grep -n "^marker_code$" fields_nl.txt | head -n1 | cut -d: -f1)
PROCESSID_ID=$(grep -n "^processid$" fields_nl.txt | head -n1 | cut -d: -f1)
TAXID_ID=$(grep -n "^taxid$" fields_nl.txt | head -n1 | cut -d: -f1)
K_ID=$(grep -n "^kingdom$" fields_nl.txt | head -n1 | cut -d: -f1)
P_ID=$(grep -n "^phylum$" fields_nl.txt | head -n1 | cut -d: -f1)
C_ID=$(grep -n "^class$" fields_nl.txt | head -n1 | cut -d: -f1)
O_ID=$(grep -n "^order$" fields_nl.txt | head -n1 | cut -d: -f1)
F_ID=$(grep -n "^family$" fields_nl.txt | head -n1 | cut -d: -f1)
G_ID=$(grep -n "^genus$" fields_nl.txt | head -n1 | cut -d: -f1)
S_ID=$(grep -n "^species$" fields_nl.txt | head -n1 | cut -d: -f1)
IDM_ID=$(grep -n "^identification_method$" fields_nl.txt | head -n1 | cut -d: -f1)
NUC_ID=$(grep -n "^nuc$" fields_nl.txt | head -n1 | cut -d: -f1)
INSDC_ACS_ID=$(grep -n "^insdc_acs$" fields_nl.txt | head -n1 | cut -d: -f1)
SRS_ID=$(grep -n "^sequence_run_site$" fields_nl.txt | head -n1 | cut -d: -f1)

rm -f out.tsv && touch out.tsv
IFS=","
while read TAXON RANK ; do 

    RANK_ID=$(grep -n "^${RANK}$" fields_nl.txt | head -n1 | cut -d: -f1)

    awk \
        -v FS="\t" \
        -v OFS="\t" \
        -v marker_id="$MARKER_CODE_ID" \
        -v processid_id="$PROCESSID_ID" \
        -v taxid_id="$TAXID_ID" \
        -v k_id="$K_ID" \
        -v p_id="$P_ID" \
        -v c_id="$C_ID" \
        -v o_id="$O_ID" \
        -v f_id="$F_ID" \
        -v g_id="$G_ID" \
        -v s_id="$S_ID" \
        -v idm_id="$IDM_ID" \
        -v nuc_id="$NUC_ID" \
        -v insdc_acs_id="$INSDC_ACS_ID" \
        -v srs_id="$SRS_ID" \
        -v rank_id="$RANK_ID" \
        -v taxon="$TAXON" \
        -v marker_code="$MARKER" \
        -v min_len="$MIN_LEN" \
        -v max_len="$MAX_LEN" \
        ' 
        # print headers
        NR == 1 { print $processid_id "\t" $taxid_id "\t" $k_id "\t" $p_id "\t" $c_id "\t" $o_id "\t" $f_id "\t" $g_id "\t" $s_id "\t" $nuc_id "\t" $idm_id "\t" $insdc_acs_id "\t" $srs_id "\t" }
        # print rows only if taxa match
        NR > 1 && $marker_id == marker_code && $rank_id == taxon { 
            # get nucleotide string without gap characters
            gsub(/-/, "", $nuc_id)
            # only print record if nucleotide string is between min and max
            if (length($nuc_id) >= min_len && length($nuc_id) <= max_len) print $processid_id "\t" $taxid_id "\t" $k_id "\t" $p_id "\t" $c_id "\t" $o_id "\t" $f_id "\t" $g_id "\t" $s_id "\t" $nuc_id "\t" $idm_id "\t" $insdc_acs_id "\t" $srs_id "\t"
        }
        ' \
        < $BOLD_TSV \
        >> out.tsv

done < taxa.csv

## remove duplicate records that may have been fetched by multiple taxonomic targets
# get header
head -n1 out.tsv > header.tsv
# get body, removing duplicate rows, and adding header back
tail -n+2 out.tsv | sort -u | cat header.tsv - > bold_extracted.tsv

# remove unneeded files
rm -f header.tsv
rm -f out.tsv
rm -f *.txt
rm -f taxa.csv