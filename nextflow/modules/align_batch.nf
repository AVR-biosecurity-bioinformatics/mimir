process ALIGN_BATCH {
    def module_name = "align_batch"
    tag "-"
    label "medium"
    container "staphb/clustalo:1.2.4"

    input:
    path(fasta_file)

    output: 
    path("*.aligned.fasta"),             emit: aligned_fasta

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