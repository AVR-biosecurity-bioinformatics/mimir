process MAKE_BLAST_DATABASE {
    def module_name = "make_blast_database"
    // tag "-"
    container "ncbi/blast:2.17.0"

    input:
    path(target_fasta)

    output: 
    path("blast_db.*"),                                         emit: blast_db

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
        ${target_fasta}
        
    """

}