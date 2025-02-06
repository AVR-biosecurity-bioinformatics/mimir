process PARSE_MARKER {
    def module_name = "parse_marker"
    tag "-"
    label "very_small"
    container "ellerbrock/alpine-bash-curl-ssl:0.3.0"

    input:
    val(marker)

    output: 
    env(GENBANK),               emit: genbank_query
    env(BOLD),                  emit: bold_query
    env(CODING),                emit: coding
    env(TYPE),                  emit: type
    path("full.hmm"),           emit: phmm
    path("seed.stockholm"),     emit: seed

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash

    ### run module code
    source ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        ${marker}
    
    """
}