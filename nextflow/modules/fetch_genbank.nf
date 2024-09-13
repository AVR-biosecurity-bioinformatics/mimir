process FETCH_GENBANK {
    def module_name = "fetch_genbank"
    tag "-"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    val(taxon)
    val(db_file)

    output: 
    tuple val(taxon), val("genbank"), path("*_genbank.rds"),                  emit: seqs
    path("*.fasta"),                                                           emit: fasta

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    taxon =                 "${taxon}"
    db_file =               "${db_file}"
    task_index =            "${task.index}"

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