process CREATE_SYNTHETIC_LARGE {
    def module_name = "create_synthetic_large"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple path(fasta_files, name: 'input*.fasta'), path(clusters_tsv)
    path(rf_counts_tsv)
    path(small_reps_txt)

    output: 
    path("synthetic_genera.fasta"),                          emit: synthetic_fasta // sequences renamed as synthetic genera, one file
    path("reps.txt"),                                        emit: reps // representative sequence names for synthetic genera
    path("counts_renamed.tsv"),                              emit: counts_renamed_tsv // counts table for renamed sequences 

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fasta_files =                    "${fasta_files}"
    clusters_tsv =                   "${clusters_tsv}"
    rf_counts_tsv =                  "${rf_counts_tsv}"
    small_reps_txt =                 "${small_reps_txt}"

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