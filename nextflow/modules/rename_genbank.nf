process RENAME_GENBANK {
    def module_name = "rename_genbank"
    tag "-"
    // label "small"
    time '5.m'
    memory '4.GB'
    cpus 1
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    val(fasta_file)
    val(ncbi_rankedlineage_noname)

    output: 
    path("renamed.fasta"),                    emit: fasta

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fasta_file =                    "${fasta_file}"
    ncbi_rankedlineage_noname =           "${ncbi_rankedlineage_noname}"

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