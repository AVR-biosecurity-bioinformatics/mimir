process SEQUENCE_TRACKER {
    def module_name = "sequence_tracker"
    tag "-"
    time '30.m'
    memory '16.GB'
    cpus 1
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(source_fates_file, name: 'sources_fates.csv')

    output: 
    path("fates_plot.pdf"),                     emit: fates_plot
    path("sf_meta.csv"),                        emit: sf_meta
    path("taxa_summary.csv"),                   emit: taxa_summary

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
      
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    sources_fates_file = "sources_fates.csv"

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