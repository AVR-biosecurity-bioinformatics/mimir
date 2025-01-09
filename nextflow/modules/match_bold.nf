process MATCH_BOLD {
    def module_name = "match_bold"
    // cache 'lenient'
    tag "-"
    // label "small"
    time '30.m'
    memory '4.GB'
    cpus 1
    container "jackscanlan/piperline-multi:0.0.1"
    fair true

    input:
    path(seq_tibble)
    path(ncbi_lineageparents)
    path(ncbi_synonyms)

    output: 
    tuple path("bold_seqs.*.fasta"), path("matching_taxids.*.csv"), path("synchanges.*.csv"),    emit: matching_data

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    seq_tibble =                    "${seq_tibble}"
    ncbi_lineageparents =           "${ncbi_lineageparents}"
    ncbi_synonyms =                 "${ncbi_synonyms}"

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
    if ("${params.rdata}" == "true") { save.image(file = "${task.process}_${task.index}.rda") } 
    })

    """
}