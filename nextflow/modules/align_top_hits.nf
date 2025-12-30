process ALIGN_TOP_HITS {
    def module_name = "align_top_hits"
    // tag "-"
    container "staphb/mafft:7.526"

    input:
    path(top_hits_tsv)
    path(seqs_fasta)

    output: 
    path("*.aligned.fasta"),                            emit: fasta

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash

    ### run module code #
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        "${top_hits_tsv}" \
        "${seqs_fasta}"
            
    """
}