process MAKE_GENCODES {
    def module_name = "make_gencodes"
    tag "-"
    // label "small"
    time '5.m'
    memory '4.GB'
    cpus 1
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(ncbi_rankedlineage)
    path(ncbi_nodes)

    output: 
    path("ncbi_gencodes.rds"),             emit: gencodes
    
    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ##  input channel variables

    ## global variables
    projectDir = "$projectDir"
    params_dict = "$params"
    ncbi_rankedlineage = "$ncbi_rankedlineage"
    ncbi_nodes = "$ncbi_nodes"

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