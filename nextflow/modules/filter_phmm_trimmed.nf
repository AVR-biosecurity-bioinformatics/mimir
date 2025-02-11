process FILTER_PHMM_TRIMMED {
    def module_name = "filter_phmm_trimmed"
    tag "-"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple path(fasta_file), path(hmmer_output, name: 'hmmer_domtblout.txt')
    val(hmm_max_evalue)
    val(hmm_min_score)
    val(hmm_max_hits)
    val(hmm_min_acc)
    val(hmm_max_gap)
    val(min_length_trimmed)

    output: 
    path("retained_trimmed.fasta"),                       emit: fasta  
    path("removed_trimmed.fasta"),                        emit: removed_fasta
    path("removed_trimmed.csv"),                          emit: removed_csv


    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ### input channel variables
    fasta_file =            "${fasta_file}"
    hmmer_output =          "${hmmer_output}"
    hmm_max_evalue =        "${hmm_max_evalue}"
    hmm_min_score =         "${hmm_min_score}"
    hmm_max_hits =          "${hmm_max_hits}"
    hmm_min_acc =           "${hmm_min_acc}"
    hmm_max_gap =           "${hmm_max_gap}"
    min_length_trimmed =    "${min_length_trimmed}"

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