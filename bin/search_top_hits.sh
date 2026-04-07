#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = query_fasta
# $4 = thresholds_csv
# $5 = n_top_hits (number of top hits to return per sequence)

THRESHOLDS_CSV=$4
N_TOP_HITS=$5

# replace each space in sequence headers of fasta files with the string "!?!?"
sed '/^>/ s/ /!?!?/g' $3 > query.fasta

echo "Starting lastal search..."
# search
lastal \
    -P 1 \
    -m 10000 \
    -k 2 \
    -u 3 \
    -D 100 \
    -N 10000 \
    -f BlastTab \
    DB \
    query.fasta \
    > results_pre.tsv
echo "Finished lastal search!"

# keep fields 1, 2 and 3, convert to original sequence names, then process to use threshold data
grep '^#' -v results_pre.tsv | \
    cut -f1,2,3 |
    sed -e 's/!?!?/ /g' |
    awk  -v RS="\n" -v ORS="\n" -v FS="\t" -v OFS="\t" '
        BEGIN {
            prev_query = ""
            file_index = 0
        }
        { 
        if ( $1 != prev_query ) {
            file_index++
        }
        print $0 > file_index".split" 
        prev_query = $1
    }' - 
echo "Finished splitting results file by query!"

## assign max LSR thresholds from threshold csv as variables in format 'max_<RANK>', in % form (ie. x100; 0.9 = 90%)
head -n1 $THRESHOLDS_CSV | tr ',' '\n' > thresholds_header.txt # header as one-column-per-line file
LSR_ID=$(grep -n "^lsr$" thresholds_header.txt | head -n1 | cut -d: -f1)
MED_ID=$(grep -n "^med$" thresholds_header.txt | head -n1 | cut -d: -f1)
MAD_ID=$(grep -n "^mad$" thresholds_header.txt | head -n1 | cut -d: -f1)
MAX_ID=$(grep -n "^max$" thresholds_header.txt | head -n1 | cut -d: -f1)
MIN_ID=$(grep -n "^min$" thresholds_header.txt | head -n1 | cut -d: -f1)

for i in family order class phylum kingdom; do
    MAX_VALUE=$( awk -v FS="," -v lsr_id="$LSR_ID" -v max_id="$MAX_ID" -v rank_name="$i" \
        'NR > 1 && $lsr_id == rank_name {print $max_id * 100}' \
        $THRESHOLDS_CSV )
    # note that the lack of " " around "=" in the next two lines is critical
    declare max_${i}="$MAX_VALUE"
    test_out=max_${i}
    # check variable is set
    if [ -z "${!test_out}" ]; then echo "max_${i} is unset" && exit 1; else echo "max_${i} = ${!test_out}"; fi
done
echo "Finished parsing thresholds!"

## work through files, keeping only the N best top hits that also have a %ID > ( relevant max threshold - 1)
echo "Starting per-hit filtering..."
touch results_filtered.tsv
for j in *.split; do 
    echo "${j}: starting..."
    ## get the query's lineage string, break into variables and determine the 0-based index of lowest classified rank (LCRI)
    # get query lineage as an array (0-6)
    IFS=$'\t' read -r -a query_taxa <<< "$(head $j -n1 | cut -f1 | sed -e 's/;/\t/g' | cut -f2,3,4,5,6,7,8)"
    echo "${j}: made query_taxa array"
    # get query lineages for each rank ('K', 'K;P', 'K;P;C' etc.)
    QL_K="${query_taxa[0]}"
    QL_P="${query_taxa[0]};${query_taxa[1]}"
    QL_C="${query_taxa[0]};${query_taxa[1]};${query_taxa[2]}"
    QL_O="${query_taxa[0]};${query_taxa[1]};${query_taxa[2]};${query_taxa[3]}"
    QL_F="${query_taxa[0]};${query_taxa[1]};${query_taxa[2]};${query_taxa[3]};${query_taxa[4]}"
    QL_G="${query_taxa[0]};${query_taxa[1]};${query_taxa[2]};${query_taxa[3]};${query_taxa[4]};${query_taxa[5]}"
    QL_S="${query_taxa[0]};${query_taxa[1]};${query_taxa[2]};${query_taxa[3]};${query_taxa[4]};${query_taxa[5]};${query_taxa[6]}"
    echo "${j}: determined query taxonomy"
    # determine highest unclassified rank index (HURI)
    for i in "${!query_taxa[@]}"; do
        if [[ "${query_taxa[$i]}" =~ "Unclassified" ]]; then
            HURI="${i}"; break
        else
            HURI="-1"
        fi
    done
    # determine lowest classified rank index (LCRI); -1 means no rank is classified
    if [[ $HURI == -1 ]]; then
        LCRI=6
    else 
        LCRI=$(( HURI - 1 ))
    fi
    echo "${j}: determined query LCR"
    # pass query info (including LCRI and the rank lineages) to awk, finding hits that pass (LSR threshold - 1)% (to allow for leniency in local vs global alignment)
    # also ignores hits that have LSR of species or genus
    awk -v RS="\n" -v ORS="\n" -v FS="\t" -v OFS="\t" \
        -v n_top_hits="$N_TOP_HITS" \
        -v max_family="$max_family" \
        -v max_order="$max_order" \
        -v max_class="$max_class" \
        -v max_phylum="$max_phylum" \
        -v max_kingdom="$max_kingdom" \
        -v q_lcri="$LCRI" \
        -v ql_k="$QL_K" \
        -v ql_p="$QL_P" \
        -v ql_c="$QL_C" \
        -v ql_o="$QL_O" \
        -v ql_f="$QL_F" \
        -v ql_g="$QL_G" \
        -v ql_s="$QL_S" \
        'BEGIN { hit_count = 1 }
        {
        # split target header into an array (1 is seqid + taxid)
        split($2, target_array, ";" )
        # determine lowest classified rank of target lineage; note that index 2 is kingdom --> index 8 is species
        for (i = 2; i <= 8; i++){
            if (target_array[i] ~ /Unclassified/)
                {t_huri = i - 2; break}
            else 
                {t_huri = -1}
        }
        if ( t_huri == -1 )
            {t_lcri = 6}
        else 
            {t_lcri = t_huri - 1}
        # determine lowest classified rank index across either sequence (lowest comparison point)
        if ( q_lcri > t_lcri )
            {lcri = t_lcri}
        else 
            {lcri = q_lcri}
        # determine lowest shared rank index
        if (ql_s == target_array[2]";"target_array[3]";"target_array[4]";"target_array[5]";"target_array[6]";"target_array[7]";"target_array[8])
            {next}
        else if (ql_g == target_array[2]";"target_array[3]";"target_array[4]";"target_array[5]";"target_array[6]";"target_array[7])
            {next}
        else if (ql_f == target_array[2]";"target_array[3]";"target_array[4]";"target_array[5]";"target_array[6])
            {lsri = 4}
        else if (ql_o == target_array[2]";"target_array[3]";"target_array[4]";"target_array[5])
            {lsri = 3}
        else if (ql_c == target_array[2]";"target_array[3]";"target_array[4])
            {lsri = 2}
        else if (ql_p == target_array[2]";"target_array[3])
            {lsri = 1}
        else if (ql_k == target_array[2])
            {lsri = 0}
        else 
            {lsri = -1 }
        # determine if lsri is lower than lcri
        if (lsri < lcri)
            # comparison is possible
            {
                # determine threshold (lsr pid - 1) for relevant rank
                if (lsri == 4)
                    {pid_threshold = (max_family - 1)}
                else if (lsri == 3)
                    {pid_threshold = (max_order - 1)}
                else if (lsri == 2)
                    {pid_threshold = (max_class - 1)}
                else if (lsri == 1)
                    {pid_threshold = (max_phylum - 1)}
                else if (lsri == 0)
                    {pid_threshold = (max_kingdom - 1)}
                else
                    {next}
                # print query and target only if threshold is exceeded
                if ($3 > pid_threshold && hit_count < n_top_hits)
                    {print $1 "\t" $2; hit_count++}
                else 
                    {next}
            }
    }' $j \
    >> results_filtered.tsv

    rm -f $j query.txt
    echo "${j}: filtered to hits above max thresholds"
done 
echo "Finished per-hit filtering!"
echo "Cleaning up..."
rm -f results_pre.tsv