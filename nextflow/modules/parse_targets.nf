process PARSE_TARGETS {
    def module_name = "parse_targets"
    tag "-"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    val(taxon)
    val(taxon_rank)
    val(entrez_key)
    val(ncbi_synonyms_file)

    output: 
    tuple path("*_name.txt"), val(taxon_rank),                      emit: taxon_name
    path("*_ncbi_id.txt"),                                          emit: ncbi_id
    path("*_bold_names.txt"),                                       emit: bold_names
    path("*_bold_ids.txt"),                                         emit: bold_ids
    path("*_bold_rank.txt"),                                        emit: bold_rank

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    taxon =                     "${taxon}"
    taxon_rank =                "${taxon_rank}"
    entrez_key =                "${entrez_key}"
    ncbi_synonyms_file =        "${ncbi_synonyms_file}"

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