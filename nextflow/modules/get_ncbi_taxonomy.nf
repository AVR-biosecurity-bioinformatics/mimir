process GET_NCBI_TAXONOMY {
    def module_name = "get_ncbi_taxonomy"
    tag "-"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:

    output: 
    path("ncbi_rankedlineage.rds"),             emit: rankedlineage
    path("ncbi_nodes.rds"),                     emit: nodes
    path("ncbi_taxidnames.rds"),                emit: taxidnames
    path("ncbi_names.rds"),                     emit: names
    path("ncbi_taxid_name_rank.rds"),           emit: taxidnamerank
    path("ncbi_synonyms.rds"),                  emit: synonyms
    path("./ncbi_taxdump"),                     emit: db_path

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables

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