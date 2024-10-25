#!/usr/bin/env nextflow
process QUERY_GENBANK {
    def module_name = "query_genbank"
    tag "-"
    cpus 1
    time 2.h
    memory 4.GB
    container "emehinovic72/edirect:latest"

    input:
    tuple val(taxon_id), val(taxon_rank)
    val(marker)

    output:
    path("*_seqs.txt"),              emit: seq_acc

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
        "${taxon_id}" \
        "${taxon_rank}" \
        "${marker}" \
        "${params.min_length}" \
        "${params.max_length}" \
        "${params.use_mito}"
    
    """

}