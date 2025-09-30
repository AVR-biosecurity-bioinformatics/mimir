#!/usr/bin/env nextflow
process QUERY_GENBANK {
    def module_name = "query_genbank"
    // tag "-"
    maxForks 10
    container "emehinovic72/edirect:latest"

    input:
    tuple val(taxon_id), val(taxon_rank)
    val(marker)
    val(min_length)
    val(max_length)
    val(use_mito)

    output:
    path("*_seqs.txt"),              emit: seq_acc

    // publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

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
        "${min_length}" \
        "${max_length}" \
        "${use_mito}"
    
    """

}