process FILTER_REDUNDANT {
    def module_name = "filter_redundant"
    // tag "-"
    container "nanozoo/bbmap:38.86--9ebcbfa"

    input:
    path(fasta_file)

    output: 
    tuple path("*.retained.fasta"), path("counts.tsv"),         emit: fasta
    path("removed.fasta"),                                      emit: removed

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash

    #### run module code
    source ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        "${fasta_file}" 
    
    """
}