process STOCKHOLM_TO_FASTA {
    def module_name = "stockholm_to_fasta"
    tag "-"
    label "very_small"
    container "library/perl:5.10.1-buster"

    input:
    path(alignment, name: 'alignment.sth')

    output: 
    path("alignment.fasta"), emit: fasta
 
    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    """
    #!/usr/bin/env bash

    # any2fasta v0.4.2 from https://github.com/tseemann/any2fasta/commit/cdec1f06c580c966a7a18874bf5c6e048bbe597f
    # in ./bin 

    ### run module code
    any2fasta.pl alignment.sth > alignment.fasta
    
    """
}