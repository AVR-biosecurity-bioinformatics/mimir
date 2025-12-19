process ALIGN_BATCH {
    def module_name = "align_batch"
    // tag "-"
    container "staphb/mafft:7.526"

    input:
    path(fasta_files)
    val(file_type)

    output: 
    path("*.aligned.fasta"),             emit: fasta

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
        "${fasta_files}" \
        ${file_type}
            
    """
}