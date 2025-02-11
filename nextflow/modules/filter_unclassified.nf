process FILTER_UNCLASSIFIED {
    def module_name = "filter_unclassified"
    tag "-"
    label "very_small"
    container "staphb/seqkit:2.8.2"

    input:
    path(fasta_file)
    val(remove_unclassified)

    output: 
    path("filter_unclassified.fasta"),      emit: fasta

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
        ${remove_unclassified}
    
    """
}