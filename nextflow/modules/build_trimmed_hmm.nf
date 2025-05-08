process BUILD_TRIMMED_HMM {
    def module_name = "build_trimmed_hmm"
    // tag "-"
    container "staphb/hmmer:3.4"

    input:
    path('seed.fasta')

    output: 
    path("trimmed.hmm"), emit: hmm 

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    """
    #!/usr/bin/env bash

    ### run module code
    hmmbuild trimmed.hmm seed.fasta
    
    """
}