process TRIM_HMM_SEED {
    def module_name = "trim_hmm_seed"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(fasta_file)
    path(primers_file)
    val(pad_fwd)
    val(pad_rev)

    output: 
    path("seed_trimmed.fasta"),                             emit: fasta
    path("primer_info.csv"),                                emit: primer_info

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fasta_file =             "${fasta_file}"
    primers_file =           "${primers_file}"
    pad_fwd =                "${pad_fwd}"
    pad_rev =                "${pad_rev}"
    
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