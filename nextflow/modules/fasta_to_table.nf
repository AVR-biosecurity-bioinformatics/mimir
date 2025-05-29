process FASTA_TO_TABLE {
    def module_name = "fasta_to_table"
    // tag "-"
    container "cicirello/gnu-on-alpine:3.20.3"

    input:
    path(fasta_list, name: '*')
    val(data_type)

    output: 
    path("output.tsv"),                            emit: tsv

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    """
    #!/usr/bin/env bash

    ### run module code

    # make header based on "data_type"
    if [[ "$data_type" == "sources" ]]; then 
        echo "name\tsource" > header.tsv
    elif [[ "$data_type" == "fates" ]]; then 
        echo "name\tfate" > header.tsv
    else 
        echo "Data type not valid"
        exit 1
    fi

    # loop through input .fasta files, appending headers and file basename into a .csv file
    touch output_noheader.tsv
    for i in $fasta_list; do
        BASENAME=\$( basename \$i .fasta )
        awk -v BASENAME="\$BASENAME" -v RS=">" -v ORS="\n" -v FS="[\r\n]+" -v OFS="\n" 'NR>1 { print \$1 "\t" BASENAME }' \$i >> output_noheader.tsv
    done

    # append header to file, then sort, ignoring header
    cat header.tsv output_noheader.tsv | \
        awk -v tb="'\t'" 'NR<2{print \$0;next}{ print \$0| " sort -k 1 -t " tb }' \
        > output.tsv

        
    """

}