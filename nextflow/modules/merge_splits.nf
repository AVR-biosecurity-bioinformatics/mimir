process MERGE_SPLITS {
    def module_name = "merge_splits"
    // tag "-"
    container "cicirello/gnu-on-alpine:3.20.3"

    input:
    path('dir*/*')

    output: 
    path("merged/*.lineage.fasta"),                            emit: fasta

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash

    ### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} 

    """

}