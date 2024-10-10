#!/usr/bin/env nextflow
process QUERY_GENBANK {
    def module_name = "query_genbank"
    tag "-"
    // label "very_high"
    container "emehinovic72/edirect:latest"

    input:
    tuple val(taxon_name), val(taxon_rank)
    val(marker)

    output:
    path("*.txt"),              emit: seq_acc

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
        "${taxon_name}" \
        "${taxon_rank}" \
        "${marker}" \
        "${params.min_length}" \
        "${params.max_length}"
    
    """

}