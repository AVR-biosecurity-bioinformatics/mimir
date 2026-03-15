process RECHECK_GENERA {
    def module_name = "recheck_genera"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(intragenus_csv)
    path(seqs_all_fasta)
    path(seqs_max_fasta)
    path(aligned_genera_cfa)
    path(synthetic_reps_txt)
    path(rf_counts_tsv)
    val(con_min_n)
    val(con_min_prop)

    output: 
    path("cg_all.csv"),                                     emit: csv
    path("cg_list.txt"),                                    emit: cg_list
    path("rf_counts_new.tsv"),                              emit: counts_tsv
    path("out.fasta"),                                      emit: fasta
    
    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
     
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    intragenus_csv =                 "${intragenus_csv}"
    seqs_all_fasta =                 "${seqs_all_fasta}"
    seqs_max_fasta =                 "${seqs_max_fasta}"
    aligned_genera_cfa =             "${aligned_genera_cfa}"
    synthetic_reps_txt =             "${synthetic_reps_txt}"
    rf_counts_tsv =                  "${rf_counts_tsv}"
    con_min_n =                      "${con_min_n}"
    con_min_prop =                   "${con_min_prop}"
    
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