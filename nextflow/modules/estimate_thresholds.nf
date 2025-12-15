process ESTIMATE_THRESHOLDS {
    def module_name = "estimate_thresholds"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(summary_csv)
    val(min_k)
    val(max_k)

    output: 
    path("thresholds.csv"),                           emit: csv
    path("estimation_plot.pdf"),                         emit: estimation_plot
    path("sample_size_histogram.pdf"),                   emit: hist

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    summary_csv =                   "${summary_csv}"
    min_k =                         "${min_k}"
    max_k =                         "${max_k}"

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