process CLUSTER_SEQUENCES {
    def module_name = "cluster_sequences"
    tag "-"
    // label "medium"
    time '2.h'
    memory '16.GB'
    // memory { 1.GB + ( 1.KB *  ) }
    cpus 4
    container "nanozoo/mmseqs2:14.7e284--11077ba"

    input:
    path(fasta_file)

    output: 
    path("clusters.tsv"),                            emit: tsv

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
        ${task.memory.getKilo()} \
        ${fasta_file} \
        ${params.cluster_threshold}
        
    """

}