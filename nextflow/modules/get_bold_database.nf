#!/usr/bin/env nextflow
process GET_BOLD_DATABASE {
    def module_name = "get_bold_database"
    tag "-"
    label "small"
    // container "cicirello/gnu-on-alpine:3.20.3"
    container "ellerbrock/alpine-bash-curl-ssl:0.3.0"

    input:
    val(bold_db_path)
    val(bold_db_url)

    output:
    tuple path("BOLD_Public*.tsv"), path("BOLD_Public.*.datapackage.json"), emit: bold_db

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'symlink'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash

    ### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        ${bold_db_path} \
        "${bold_db_url}"
    
    """

}