process TRIM_TO_PRIMERS {
    def module_name = "trim_to_primers"
    tag "-"
    // label "medium"
    time '1.h'
    memory '8.GB'
    cpus 1
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(fasta_file)
    val(primer_fwd)
    val(primer_rev)
    val(remove_primers)
    val(max_primer_mismatches)
    val(min_length_trimmed)

    output: 
    path("trimmed.fasta"),                                  emit: fasta
    path("trimmed.removed.fasta"),                          emit: removed

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fasta_file =                    "${fasta_file}"
    primer_fwd =                    "${primer_fwd}"
    primer_rev =                    "${primer_rev}"
    remove_primers =                "${remove_primers}"
    max_primer_mismatches =         "${max_primer_mismatches}"
    min_length_trimmed =            "${min_length_trimmed}"

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