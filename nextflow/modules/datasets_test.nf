#!/usr/bin/env nextflow
process DATASETS_TEST {
    def module_name = "datasets_test"
    tag "-"
    // label "very_high"
    container "staphb/ncbi-datasets:16.22.1"

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