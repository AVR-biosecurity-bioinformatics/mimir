#!/usr/bin/env nextflow
process TEMPLATE_BASH {
    def module_name = "template_bash"
    tag "-"
    // label "very_high"
    // container "staphb/ncbi-datasets:16.22.1"

    input:

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
        ${task.cpus} 
    
    """

}