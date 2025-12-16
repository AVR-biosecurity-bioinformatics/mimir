process SPLIT_BY_CLASSIFICATION {
    def module_name = "split_by_classification"
    // tag "-"
    container "cicirello/gnu-on-alpine:3.20.3"

    input:
    path(fasta_files)

    output: 
    path("*.full.fasta"),                            emit: full
    path("*.part.fasta"),                            emit: part


    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash

    ### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        "${fasta_files}"
        
    """

}