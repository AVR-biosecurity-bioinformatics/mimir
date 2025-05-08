process JOIN_SOURCES_FATES {
    def module_name = "join_sources_fates"
    // tag "-"
    container "cicirello/gnu-on-alpine:3.20.3"

    input:
    path('sources.csv')
    path('fates.csv')

    output: 
    path("sources_fates.csv"),                            emit: csv

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    """
    #!/usr/bin/env bash

    ### run module code
    join sources.csv fates.csv -t , --header > sources_fates.csv
        
    """

}