process ALIGN_PRIMERS_TO_SEED {
    def module_name = "align_primers_to_seed"
    tag "-"
    time '30.m'
    memory '8.GB'
    cpus 1
    container "staphb/mafft:7.526"

    input:
    path(fasta_file)
    path(primers)

    output: 
    path("aligned_primers.fasta"),             emit: fasta

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
        "${fasta_file}" \
        ${primers}
    
    """
}