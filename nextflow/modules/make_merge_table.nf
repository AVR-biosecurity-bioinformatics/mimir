process MAKE_MERGE_TABLE {
    def module_name = "make_merge_table"
    tag "-"
    label "medium"
    // container "library/ruby:3.3.7-alpine3.21"
    container "nanozoo/ruby:3.2.2--fdadcb3"

    input:
    path(fasta_files, name: 'alignment*.fasta')

    output: 
    tuple path("combined.fasta"),  path("subMSAtable"),           emit: fasta_table

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    // from: https://mafft.cbrc.jp/alignment/software/makemergetable.rb
    """
    #!/usr/bin/env bash

    ### make merge table
    makemergetable.rb \
        $fasta_files \
        > subMSAtable
    
    # throw error if subMSAtable is empty
    if [ -s subMSAtable ]; then
        ### concatenate files together
        cat $fasta_files > combined.fasta
    else 
        echo "subMSAtable is empty"
        exit 1
    fi

    """
}