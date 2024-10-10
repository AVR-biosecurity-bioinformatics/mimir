process RENAME_GENBANK {
    def module_name = "rename_genbank"
    tag "-"
    time 5.m
    cpus 1
    memory 1.GB
    container "staphb/seqkit:2.8.2"

    input:
    tuple path(fasta), path(taxids)

    output: 
    path("*renamed.fasta"),                            emit: fasta

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
        ${fasta} \
        ${taxids}
        
    """

}