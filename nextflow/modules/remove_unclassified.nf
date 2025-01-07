process REMOVE_UNCLASSIFIED {
    def module_name = "remove_unclassified"
    tag "${task.index}"
    label "very_small"
    container "staphb/seqkit:2.8.2"

    input:
    path(fasta_file)
    val(remove_unclassified)

    output: 
    path("remove_unclassified.*.fasta"),      emit: fasta

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash

    ### run module code
    source ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        "${fasta_file}" \
        ${remove_unclassified} \
        ${task.index}
    
    """
}