process JOIN_SOURCES_FATES {
    def module_name = "join_sources_fates"
    // tag "-"
    container "cicirello/gnu-on-alpine:3.20.3"

    input:
    path('sources.tsv')
    path('fates.tsv')

    output: 
    path("sources_fates.tsv"),                            emit: tsv

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    """
    #!/usr/bin/env bash

    ### run module code
    join sources.tsv fates.tsv -t \$'\t' --header > sources_fates.tsv
        
    """

}