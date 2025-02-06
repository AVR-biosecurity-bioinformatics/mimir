process ALIGN_CORE {
    def module_name = "align_core"
    tag "-"
    // label "very_high"
    time '4.h'
    memory '16.GB'
    cpus 16
    container "staphb/mafft:7.526"


    input:
    path('core.fasta')

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
        core.fasta
    
    """
}