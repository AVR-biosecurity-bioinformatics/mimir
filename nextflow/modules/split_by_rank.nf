process SPLIT_BY_RANK {
    def module_name = "split_by_rank"
    tag "-"
    // label "medium"
    // time '1.h'
    time { 10.m + ((fasta_file.size() / 100000) * 1.s ) }
    // 10 minutes + 1 second for every 100KB file size
    memory '2.GB'
    // memory { 1.GB + ( 1.KB *  ) }
    cpus 1
    container "cicirello/gnu-on-alpine:3.20.3"

    input:
    path(fasta_file)
    val(split_rank)

    output: 
    path("*.lineage.fasta"),                            emit: fasta

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    // println ( task.time )
    """
    #!/usr/bin/env bash

    ### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        ${fasta_file} \
        ${split_rank} \
        ${task.time.toSeconds()}
        
    """

}