process CONCAT_CSV {
    def module_name = "concat_csv"
    tag "-"
    // label "medium"
    time '30.m'
    memory '2.GB'
    cpus 1
    container "cicirello/gnu-on-alpine:3.20.3"

    input:
    path(fasta_list, name: '*')
    val(data_type)

    output: 
    path("output.csv"),                            emit: csv

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    """
    #!/usr/bin/env bash

    ### run module code

    # make header based on "data_type"
    if [[ "$data_type" == "sources" ]]; then 
        echo "name,source" > header.csv
    elif [[ "$data_type" == "fates" ]]; then 
        echo "name,fate" > header.csv
    else 
        echo "Data type not valid"
        exit 1
    fi

    # loop through input .fasta files, appending headers and file basename into a .csv file
    touch output.csv
    for i in $fasta_list; do
        BASENAME=\$( basename \$i .fasta )
        awk -v BASENAME="\$BASENAME" -v RS=">" -v ORS="\n" -v FS="[\r\n]+" -v OFS="\n" 'NR>1 { print \$1 "," BASENAME }' \$i >> output_noheader.csv
    done

    # append header to file, then sort, ignoring header
    cat header.csv output_noheader.csv | \
        awk 'NR<2{print \$0;next}{print \$0| "sort -k 1 -t ,"}' \
        > output.csv

        
    """

}