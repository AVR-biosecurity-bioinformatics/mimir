process FILTER_DUPLICATES {
    def module_name = "filter_duplicates"
    tag "-"
    label "very_small"
    container "staphb/seqkit:2.8.2"

    input:
    path(fasta_file)

    output: 
    path("seqs_deduplicated.fasta"),                            emit: fasta

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash

    ### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        ${fasta_file} 
        
    """

}