process FORMAT_OUTPUT {
    def module_name = "format_output"
    tag "-"
    label "medium"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    val(fasta_file)
    val(add_root)
    val(aligned_output)
    val(compressed_output)

    output: 
    path("final_database.fasta"),                           emit: fasta

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fasta_file =                   "${fasta_file}"

    ## global variables
    projectDir = "$projectDir"
    params_dict = "$params"
    add_root = "${add_root}"
    aligned_output = "${aligned_output}"
    compressed_output = "${compressed_output}"

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