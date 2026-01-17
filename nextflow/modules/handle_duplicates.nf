process HANDLE_DUPLICATES {
    def module_name = "handle_duplicates"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(bold_seqs_file, name: 'bold_in.fasta')
    path(genbank_seqs_file, name: 'genbank_in.fasta')
    path(gb_dups_file, name: 'dups_in.csv')
    path(rankedlineage_noname)
    val(prefer_source)
    val(duplicate_taxonomy)
    val(duplicate_sequence)

    output: 
    path("bold_out.fasta"),                         emit: bold
    path("genbank_out.fasta"),                      emit: genbank
    path("bold_gb_duplicates.csv"),                 emit: gb_dups

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
     
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    bold_seqs_file =                    "${bold_seqs_file}"
    genbank_seqs_file =                 "${genbank_seqs_file}"
    gb_dups_file =                      "${gb_dups_file}"
    rankedlineage_noname =              "${rankedlineage_noname}"
    prefer_source =                     "${prefer_source}"
    duplicate_taxonomy =                "${duplicate_taxonomy}"
    duplicate_sequence =                "${duplicate_sequence}"

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