process TRIM_PHMM {
    def module_name = "trim_phmm"
    tag "-"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    val(phmm_model_file)
    val(primer_fwd)
    val(primer_rev)
    val(remove_primers)

    output: 
    path("phmm_model_trimmed.rds"),                  emit: trimmed_model

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    phmm_model_file =        "${phmm_model_file}"
    primer_fwd =             "${primer_fwd}"
    primer_rev =             "${primer_rev}"
    remove_primers =         "${remove_primers}"

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
    if ("${params.rdata}" == "true") { save.image(file = "${task.process}_${task.index}.rda") } 
    })

    """
}