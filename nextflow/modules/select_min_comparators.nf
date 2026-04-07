process SELECT_MIN_COMPARATORS {
    def module_name = "select_min_comparators"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(cg_csv)
    path(cg_list)
    path(seqs_max_fasta)
    val(con_min_n)
    val(con_min_prop)

    output: 
    path("min_comparators.cfa"),                            emit: cfa
    
    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
     
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    cg_csv =                         "${cg_csv}"
    cg_list =                        "${cg_list}"
    seqs_max_fasta =                 "${seqs_max_fasta}"
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