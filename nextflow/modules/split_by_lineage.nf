process SPLIT_BY_LINEAGE {
    def module_name = "split_by_lineage"
    tag "-"
    // label "medium"
    time '1.h'
    memory '8.GB'
    cpus 1
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    val(fasta_file)
    val(split_rank_input)

    output: 
    path("*.grouped.fasta"),                                  emit: split_fasta

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fasta_file =                   "${fasta_file}"
    split_rank_input =                   "${split_rank_input}"

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