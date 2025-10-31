process FILTER_PHMM_FULL {
    def module_name = "filter_phmm_full"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple path(fasta_file), path(translations), path(hmmer_output, name: 'hmmer_domtblout.txt')
    val(phmm_max_evalue)
    val(phmm_min_score)
    val(phmm_max_hits)
    val(phmm_min_acc)
    val(phmm_max_gap)
    val(phmm_min_length)
    val(phmm_min_cov)

    output: 
    tuple path("retained.fasta"), path("translations_retained.fasta"),  emit: retained 
    path("removed.fasta"),                                              emit: removed_fasta
    path("removed.csv"),                                                emit: removed_csv


    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ### input channel variables
    fasta_file =            "${fasta_file}"
    translations =          "${translations}" 
    hmmer_output =          "${hmmer_output}"
    phmm_max_evalue =       "${phmm_max_evalue}"
    phmm_min_score =        "${phmm_min_score}"
    phmm_max_hits =         "${phmm_max_hits}"
    phmm_min_acc =          "${phmm_min_acc}"
    phmm_max_gap =          "${phmm_max_gap}"
    phmm_min_length =       "${phmm_min_length}"
    phmm_min_cov =          "${phmm_min_cov}"
       
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