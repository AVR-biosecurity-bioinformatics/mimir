process GET_MIN_COMPONENTS {
    def module_name = "get_min_components"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(comparators_cfa)
    path(cg_csv)
    path(seqs_file)
    path(counts_file)
    path(thresholds_csv)
    val(component_group_size)

    output: 
    path("component_group*.fasta"),                                  emit: fasta

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    comparators_cfa =                      "${comparators_cfa}"
    cg_csv =                               "${cg_csv}"
    seqs_file =                            "${seqs_file}"
    counts_file =                          "${counts_file}"
    thresholds_csv =                       "${thresholds_csv}"
    component_group_size =                 "${component_group_size}"
    
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