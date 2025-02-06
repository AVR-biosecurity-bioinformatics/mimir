process MERGE_ALIGNMENTS {
    def module_name = "merge_alignments"
    tag "-"
    // label "very_high"
    cpus 8
    time '2.h'
    memory '16.GB'
    // container "staphb/clustalo:1.2.4"
    container "staphb/mafft:7.526"


    input:
    tuple path('combined.fasta'), path('subMSAtable')

    output: 
    path("aligned.fasta"),             emit: aligned_fasta

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash

    ### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        combined.fasta \
        subMSAtable
    
    """
}