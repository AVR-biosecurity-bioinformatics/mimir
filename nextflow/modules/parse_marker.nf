process PARSE_MARKER {
    def module_name = "parse_marker"
    tag "-"
    label "very_small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    val(marker)

    output: 
    env(GENBANK),               emit: genbank_query
    env(BOLD),                  emit: bold_query
    env(CODING),                emit: coding

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