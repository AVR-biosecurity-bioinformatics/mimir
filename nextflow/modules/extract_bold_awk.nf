process EXTRACT_BOLD_AWK {
    def module_name = "extract_bold_awk"
    // tag "-"
    container "cicirello/gnu-on-alpine:3.20.3"

    input:
    tuple path(db_tsv_file), path(db_meta_file)
    path(bold_targets)
    val(marker)
    val(min_length_input)
    val(max_length_input)

    output: 
    path("bold_extracted.tsv"),                     emit: tsv

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
        "${db_tsv_file}" \
        "${db_meta_file}" \
        "${bold_targets}" \
        "${marker}" \
        "${min_length_input}" \
        "${max_length_input}" 
        
    """
}