process SORT_BY_LINEAGE {
    def module_name = "sort_by_lineage"
    tag "-"
    // label "medium"
    // time '1.h'
    time { 1.m + ((fasta_file.size() / 100000) * 1.s ) }
    // 1 minutes + 1 second for every 100KB file size
    memory '2.GB'
    // memory { 1.GB + ( 1.KB *  ) }
    cpus 1
    container "cicirello/gnu-on-alpine:3.20.3"
 
    input:
    path(fasta_file)

    output: 
    path("sorted.fasta"),                            emit: fasta

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
        ${task.time.toSeconds()}
        
    """

}