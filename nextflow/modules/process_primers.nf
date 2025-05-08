process PROCESS_PRIMERS {
    def module_name = "process_primers"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    val(primer_fwd)
    val(primer_rev)
    path(gencodes_file)
    val(marker_type)

    output: 
    path("primers.fasta"),                                  emit: nuc_fasta
    path("translated.fasta"),                               emit: translated

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    primer_fwd =                    "${primer_fwd}"
    primer_rev =                    "${primer_rev}"
    gencodes_file =                 "${gencodes_file}"
    marker_type =                   "${marker_type}"

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