process ALIGN_SINGLE {
    def module_name = "align_single"
    tag "-"
    label "very_high"
    container "staphb/clustalo:1.2.4"

    input:
    val(fasta_file)

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
        "${fasta_file}"
    
    """
}