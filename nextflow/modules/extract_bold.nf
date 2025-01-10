process EXTRACT_BOLD {
    def module_name = "extract_bold"
    // cache 'lenient'
    tag "-"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"
    // fair true

    input:
    tuple path(db_tsv_file), path(db_meta_file)
    path(bold_names_file)
    path(bold_rank_file)
    val(marker)

    output: 
    path("bold_db_targets.rds"),                  emit: tibble, optional: true

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    db_tsv_file =                 "${db_tsv_file}"
    db_meta_file =                "${db_meta_file}"
    bold_names_file =             "${bold_names_file}"
    bold_rank_file =              "${bold_rank_file}"
    marker =                      "${marker}"

    ## global variables
    projectDir = "$projectDir"
    params_dict = "$params"

    tryCatch({
    ### source functions and themes, load packages, and import Nextflow params
    ### from "bin/process_start.R"
    sys.source("${projectDir}/bin/process_start.R", envir = .GlobalEnv)

    ### run module code
    sys.source(
        "${projectDir}/bin/$module_script", # run script
        envir = .GlobalEnv # this allows import of existing objects like projectDir
    )
    }, finally = {
    ### save R environment for debugging
    if ("${params.rdata}" == "true") { save.image(file = "${task.process}.rda") } 
    })

    """
}