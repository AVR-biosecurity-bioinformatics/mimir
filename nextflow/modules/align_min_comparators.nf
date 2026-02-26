process ALIGN_MIN_COMPARATORS {
    def module_name = "align_min_comparators"
    // tag "-"
    container "staphb/mafft:7.526"

    input:
    path(cfa_file)

    output: 
    path("alignment.cfa"),             emit: cfa

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run module code #
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        "${cfa_file}" 
            
    """
}