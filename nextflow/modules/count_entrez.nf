#!/usr/bin/env nextflow
process COUNT_ENTREZ {
    def module_name = "count_entrez"
    tag "-"
    // label "very_high"
    container "emehinovic72/edirect:latest"

    input:
    val(gene)
    val(taxon)

    output:


    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/bash

    ### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        ${gene} \
        ${taxon}
    
    """

}