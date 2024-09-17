process REMOVE_EXACT_DUPLICATES {
    def module_name = "remove_exact_duplicates"
    tag "-"
    label "very_small"
    container "staphb/seqkit:2.8.2"

    input:
    val(fasta_file)

    output: 
    // path("seqs_combined.rds"),                  emit: seqs
    path("seqs_deduplicated.fasta"),                            emit: fasta

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/bash

    ### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        ${fasta_file} 
        
    """

}