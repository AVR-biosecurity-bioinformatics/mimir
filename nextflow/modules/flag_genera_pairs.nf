process FLAG_GENERA_PAIRS {
    def module_name = "flag_genera_pairs"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(top_hits_file)
    path(thresholds_file)
    path(seqs_file)
    path(counts_file)

    output: 
    path("flagged_genera.csv"),                                  emit: csv

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    top_hits_file =                     "${top_hits_file}"
    thresholds_file =                   "${thresholds_file}"
    seqs_file =                         "${seqs_file}"
    counts_file =                       "${counts_file}"

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