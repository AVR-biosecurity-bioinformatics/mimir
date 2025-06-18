process CLUSTER_SEQUENCES {
    def module_name = "cluster_sequences"
    // tag "-"
    container "nanozoo/mmseqs2:14.7e284--11077ba"

    input:
    path(fasta_file)
    val(cluster_threshold)

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
        ${cluster_threshold}
        
    """

}