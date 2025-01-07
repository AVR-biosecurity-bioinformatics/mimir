process RESOLVE_SYNONYMS {
    def module_name = "resolve_synonyms"
    tag "-"
    cpus 1
    memory 5.GB
    time 10.m
    // label "medium"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(seqs_file), val(seq_source)
    val(db_path)

    output: 
    tuple path("*_seqs_resolved.rds"), val(seq_source),         emit: seqs
    path("*.fasta"),                                            emit: fasta, optional: true

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    seqs_file =              "${seqs_file}"
    db_path =                "${db_path}"

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