process COMBINE_CHUNKS {
    def module_name = "combine_chunks"
    // tag "-"
    container "staphb/seqkit:2.8.2"

    input:
    path('seq???????.fasta')
    val(dealign)

    output: 
    path("chunks_combined.fasta"),                            emit: fasta

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
        'seq*.fasta' \
        ${dealign}
        
    """

}