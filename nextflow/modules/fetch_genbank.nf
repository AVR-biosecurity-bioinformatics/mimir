#!/usr/bin/env nextflow
process FETCH_GENBANK {
    def module_name = "fetch_genbank"
    tag "$acc_list"
    time 10.m
    cpus 1
    memory 1.GB
    container "emehinovic72/edirect:latest"
    maxForks 10

    input:
    path(acc_list)
    val(entrez_key)

    output:
    tuple path("*.fasta"), path("*.taxids.txt"),                  emit: fetched_seqs

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash

    ### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        ${acc_list} \
        ${entrez_key}
    
    """

}