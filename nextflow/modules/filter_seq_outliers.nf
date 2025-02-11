process FILTER_SEQ_OUTLIERS {
    def module_name = "filter_seq_outliers"
    // tag "-"
    // label "medium"
    time '1.h'
    memory '4.GB'
    cpus 1
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(fasta_files)
    val(dist_threshold)

    output: 
    path("*.retained.fasta"),                                  emit: retained_fasta
    path("*.removed.fasta"),                                   emit: removed

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fasta_files =                   "${fasta_files}"
    dist_threshold =                "${dist_threshold}"

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