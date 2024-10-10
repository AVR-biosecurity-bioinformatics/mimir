process FETCH_GENBANK_OLD {
    def module_name = "fetch_genbank_old"
    tag "$taxon"
    // label "medium"
    cpus 1
    memory 2.GB
    time {
        int seq_count = count_file.getBaseName() as int
        if ( params.entrez_key ) {
            extraTime = Duration.of(seq_count * 105) // add 105 milliseconds per sequence
        } else {
            extraTime = Duration.of(seq_count * 340) // add 340 milliseconds per sequence
        }
        5.m + extraTime
    }
    container "jackscanlan/piperline-multi:0.0.1"
    maxForks 10

    input:
    tuple val(taxon), path(count_file)
    val(db_file)

    output: 
    path("*_genbank.rds"),                  emit: seqs
    path("*.fasta"),                        emit: fasta

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