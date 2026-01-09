process BLAST_TOP_HITS {
    def module_name = "blast_top_hits"
    // tag "-"
    container "ncbi/blast:2.17.0"

    input:
    path(query_fasta)
    path(blast_db)
    val(n_top_hits)

    output: 
    path("results_filtered.tsv"),                                         emit: tsv

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
        "${query_fasta}" \
        ${n_top_hits}
        
    """

}