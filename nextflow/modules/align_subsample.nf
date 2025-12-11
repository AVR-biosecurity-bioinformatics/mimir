process ALIGN_SUBSAMPLE {
    def module_name = "align_subsample"
    // tag "-"
    container "staphb/mafft:7.526"


    input:
    path('input.fasta')

    output: 
    path("aligned.fasta"),             emit: fasta

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash

    #### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        input.fasta
    
    """
}