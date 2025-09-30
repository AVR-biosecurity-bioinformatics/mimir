#!/usr/bin/env nextflow
process FETCH_GENBANK {
    def module_name = "fetch_genbank"
    // tag "-"
    container "emehinovic72/edirect:latest"
    maxForks 10

    input:
    path(acc_list, name: 'accessions.txt')
    val(entrez_key)

    output:
    tuple path("sequences.gb"), path("accessions.txt"),                  emit: fetched_seqs

    // publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash
    
    #### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        accessions.txt \
        ${entrez_key}
    
    """

}