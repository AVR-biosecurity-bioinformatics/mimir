process FILTER_TAX_OUTLIERS {
    def module_name = "filter_tax_outliers"
    tag "-"
    // label "medium"
    time '1.h'
    memory '16.GB'
    cpus 1
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(fasta_file)
    path(cluster_tsv)
    val(db_file)

    output: 
    path("seqs_decontaminated.fasta"),                                  emit: fasta
    path("removed.fasta"),                                              emit: removed

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fasta_file =                   "${fasta_file}"
    cluster_tsv =                   "${cluster_tsv}"
    db_file =                     "${db_file}"

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