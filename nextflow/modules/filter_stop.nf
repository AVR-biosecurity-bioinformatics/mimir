process FILTER_STOP {
    def module_name = "filter_stop"
    tag "-"
    label "very_small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple path(fasta_file), val(coding), val(marker_type), path(ncbi_gencodes)

    output: 
    path("filter_stop.fasta"),                emit: fasta

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fasta_file =             "${fasta_file}"
    coding =                 "${coding}"
    marker_type =            "${marker_type}"
    ncbi_gencodes =          "${ncbi_gencodes}"

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