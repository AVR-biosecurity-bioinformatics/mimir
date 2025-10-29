process FILTER_PHMM_AMPLICON {
    def module_name = "filter_phmm_amplicon"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple path(fasta_file), path(hmmer_output, name: 'hmmer_domtblout.txt')
    path(primer_info_file)
    val(trim_to_amplicon)
    val(amplicon_min_length)
    val(amplicon_min_cov)
    val(remove_primers)

    output: 
    path("retained_amplicon.fasta"),                       emit: fasta  
    path("removed_amplicon.fasta"),                        emit: removed_fasta
    path("removed_amplicon.csv"),                          emit: removed_csv


    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
       
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fasta_file =            "${fasta_file}"
    primer_info_file =      "${primer_info_file}"
    hmmer_output =          "${hmmer_output}"
    trim_to_amplicon =      "${trim_to_amplicon}"
    amplicon_min_length =   "${amplicon_min_length}"
    amplicon_min_cov =      "${amplicon_min_cov}"
    remove_primers =        "${remove_primers}"
 
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