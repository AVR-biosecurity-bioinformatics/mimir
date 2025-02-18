process CONCAT_CSV {
    def module_name = "concat_csv"
    tag "-"
    // label "medium"
    time '30.m'
    memory '2.GB'
    cpus 1
    container "cicirello/gnu-on-alpine:3.20.3"

    input:
    path('input*.csv')

    output: 
    path("output.csv"),                            emit: csv

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    """
    #!/usr/bin/env bash

    ### run module code
    # combine files
    awk 'FNR==1 && NR!=1{next;}{print}' input*.csv > output_unsorted.csv

    # sort output file (ignoring header)
    cat output_unsorted.csv | awk 'NR<2{print \$0;next}{print \$0| "sort -k 1 -t ,"}' > output.csv
        
    """

}