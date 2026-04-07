process ALIGN_CORE {
    def module_name = "align_core"
    // tag "-"
    container "staphb/mafft:7.526"

    input:
    tuple path('core.fasta'), path(other_fasta)

    output: 
    tuple path("core_aligned.fasta"), path(other_fasta),             emit: fasta

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
        core.fasta
    
    """
}