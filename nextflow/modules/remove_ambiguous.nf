process REMOVE_AMBIGUOUS {
    def module_name = "remove_ambiguous"
    tag "-"
    label "very_small"
    container "cicirello/gnu-on-alpine:3.20.3"

    input:
    path(fasta_file)

    output: 
    path("remove_ambiguous.fasta"),      emit: fasta
    path("removed.fasta"),          emit: removed

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash

    ### run module code
    source ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        "${fasta_file}" 
    
    """
}