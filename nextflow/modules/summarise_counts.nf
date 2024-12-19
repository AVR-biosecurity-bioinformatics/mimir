process SUMMARISE_COUNTS {
    def module_name = "summarise_counts"
    tag "-"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    val(count_genbank)
    val(count_bold)
    val(count_mito)
    val(count_genome)
    val(count_internal)
    val(count_external)
    val(count_input)
    val(count_remove_unclassified)
    val(count_filter_phmm)
    val(count_filter_stop)
    val(count_remove_exact)
    val(count_remove_tax_outliers)
    val(count_remove_seq_outliers)
    val(count_prune_groups)

    output: 
    path("*.csv"),                  emit: csv
    path("*.pdf"),                  emit: pdf

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    count_genbank = "${count_genbank}"
    count_bold = "${count_bold}"
    count_mito = "${count_mito}"
    count_genome = "${count_genome}"
    count_internal = "${count_internal}"
    count_external = "${count_external}"
    count_input = "${count_input}"
    count_remove_unclassified = "${count_remove_unclassified}"
    count_filter_phmm = "${count_filter_phmm}"
    count_filter_stop = "${count_filter_stop}"
    count_remove_exact = "${count_remove_exact}"
    count_remove_tax_outliers = "${count_remove_tax_outliers}"
    count_remove_seq_outliers = "${count_remove_seq_outliers}"
    count_prune_groups = "${count_prune_groups}"

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
    if ("${params.rdata}" == "true") { save.image(file = "${task.process}_${task.index}.rda") } 
    })

    """
}