process CHUNK_TAXON {
    def module_name = "chunk_taxon"
    tag "-"
    label "long"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    val(taxon)
    val(entrez_key)
    val(chunk_rank)

    output: 
    path("tax_list.txt"),                  emit: tax_list

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    taxon =                 "${taxon}"
    entrez_key =               "${entrez_key}"
    chunk_rank =            "${chunk_rank}"

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