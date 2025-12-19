process CLUSTER_MMSEQS {
    def module_name = "cluster_mmseqs"
    // tag "-"
    container "nanozoo/mmseqs2:14.7e284--11077ba"

    input:
    path(fasta_files)
    path(thresholds_csv)
    val(process_type)

    output: 
    tuple path(fasta_files), path("clusters.tsv"),                            emit: clusters

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
        ${task.memory.getKilo()} \
        "${fasta_files}" \
        ${thresholds_csv} \
        ${process_type}
        
    """

}