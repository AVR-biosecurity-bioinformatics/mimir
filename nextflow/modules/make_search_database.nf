process MAKE_SEARCH_DATABASE {
    def module_name = "make_search_database"
    // tag "-"
    container "quay.io/biocontainers/last:1651--h7f5d12c_0"

    input:
    path(target_fasta)

    output: 
    path("DB.*"),                                         emit: db

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash
       
    ### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        ${target_fasta}
        
    """

}