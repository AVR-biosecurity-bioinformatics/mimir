process ALIGN_OTHER {
    def module_name = "align_other"
    tag "-"
    // label "very_high"
    time '8.h'
    memory '16.GB'
    cpus 16
    container "staphb/mafft:7.526"


    input:
    path('core.aligned.fasta')
    path('other.fasta')

    output: 
    path("all.aligned.fasta"),             emit: fasta

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
        core.aligned.fasta \
        other.fasta

    
    """
}