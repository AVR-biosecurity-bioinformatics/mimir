process HMMSEARCH_TRIMMED {
    def module_name = "hmmsearch_trimmed"
    // tag "-"
    container "staphb/hmmer:3.4"

    input:
    tuple path('retained_full.fasta'), path('translations.fasta')
    path('trimmed.hmm')

    output: 
    tuple path('retained_full.fasta'), path("hmmer_domtblout.txt"),        emit: hmmer_output

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
        trimmed.hmm
    
    """
}