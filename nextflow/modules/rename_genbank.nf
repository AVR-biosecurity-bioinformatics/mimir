process RENAME_GENBANK {
    def module_name = "rename_genbank"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple path(gb_file), path(accessions_file)
    path(ncbi_rankedlineage_noname)
    val(placeholder_as_unclassified)
    val(digits_as_unclassified)

    output: 
    path("renamed.fasta"),                    emit: fasta
    path("accessions_failed.txt"),            emit: accessions_failed

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
     
    ### defining Nextflow environment variables as R variables
    ### input channel variables
    gb_file =                           "${gb_file}"
    accessions_file =                   "${accessions_file}"
    ncbi_rankedlineage_noname =         "${ncbi_rankedlineage_noname}"
    placeholder_as_unclassified =       "${placeholder_as_unclassified}"
    digits_as_unclassified =            "${digits_as_unclassified}"

    ### global variables
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