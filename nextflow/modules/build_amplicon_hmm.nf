process BUILD_AMPLICON_HMM {
    def module_name = "build_amplicon_hmm"
    // tag "-"
    container "staphb/hmmer:3.4"

    input:
    path('seed.fasta')

    output: 
    path("amplicon.hmm"), emit: hmm 

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    """
    #!/usr/bin/env bash

    ### run module code
    hmmbuild amplicon.hmm seed.fasta
    
    """
}