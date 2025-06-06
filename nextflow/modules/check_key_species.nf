process CHECK_KEY_SPECIES {
    def module_name = "check_key_species"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(sf_meta_file, name: 'sf_meta.tsv')
    path(key_species_list, name: 'key_species_list.txt')

    output: 
    path("*.seqnames.txt"),                emit: seq_names
    path("key_sequences.tsv"),             emit: key_sequences_tsv
    path("key_species_summary.tsv"),       emit: summary

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
      
    ### defining Nextflow environment variables as R variables
    ### input channel variables
    sf_meta_file = "sf_meta.tsv"
    key_species_list = "key_species_list.txt"

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