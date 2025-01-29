process ALIGN_PRIMERS {
    def module_name = "align_primers"
    tag "-"
    time '1.h'
    memory '16.GB'
    cpus 8
    container "staphb/mafft:7.526"

    input:
    path(fasta_file)
    path(primers)

    output: 
    path("aligned_primers.fasta"),             emit: fasta

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/env bash
    ####
    ### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        "${fasta_file}" \
        ${primers}
    
    """
}