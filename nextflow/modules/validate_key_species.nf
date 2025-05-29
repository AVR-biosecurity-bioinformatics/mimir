process VALIDATE_KEY_SPECIES {
    def module_name = "validate_key_species"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(seq_names_file, name: 'seq_names.txt')
    path(final_db_file, name: 'final_db.fasta')
    val(add_root)

    output: 
    path("*.identical_seqs.tsv"),              emit: identical_seqs
    path("*.key_pid.tsv"),                     emit: key_pid

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
      
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    seq_names_file =        "seq_names.txt"
    final_db_file =         "final_db.fasta"
    add_root =              "${add_root}"

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