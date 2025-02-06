process HMMSEARCH_FULL {
    def module_name = "hmmsearch_full"
    tag "-"
    // label "very_high"
    time '10.m'
    memory '8.GB'
    cpus 1
    container "staphb/hmmer:3.4"

    input:
    tuple path(fasta_file), path(translations, name: 'translations.fasta')
    path(hmm, name: 'full.hmm')

    output: 
    tuple path(fasta_file), path(translations), path("hmmer_domtblout.txt"),        emit: hmmer_output

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash

    ### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        translations.fasta \
        full.hmm
    
    """
}