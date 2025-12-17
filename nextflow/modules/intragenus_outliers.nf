process INTRAGENUS_OUTLIERS {
    def module_name = "intragenus_outliers"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(fasta_files)
    path(rf_counts_tsv)
    path(thresholds_csv)
    val(con_min_n)
    val(con_min_prop)

    output: 
    path("gs_tibble.csv"),                                  emit: csv
    path("retained.fasta"),                                 emit: retained
    path("gminor.fasta"),                                   emit: removed_g
    path("sminor.fasta"),                                   emit: removed_s
    

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fasta_files =                    "${fasta_files}"
    rf_counts_tsv =                  "${rf_counts_tsv}"
    thresholds_csv =                 "${thresholds_csv}"
    con_min_n =                      "${con_min_n}"
    con_min_prop =                   "${con_min_prop}"

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