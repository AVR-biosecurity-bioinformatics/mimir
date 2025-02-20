process ALIGN_BATCH {
    def module_name = "align_batch"
    tag "-"
    label "align"
    // time '1.h'
    // memory '4.GB'
    // cpus 4
    container "staphb/clustalo:1.2.4"

    input:
    tuple val(fasta_file), path(counts_file)

    output: 
    tuple path("*.aligned.fasta"), path(counts_file),             emit: aligned_fasta

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash

    ### run module code #
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        "${fasta_file}" 
            
    """
}