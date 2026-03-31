process SEARCH_TOP_HITS {
    def module_name = "search_top_hits"
    // tag "-"
    container "quay.io/biocontainers/last:1651--h7f5d12c_0"

    input:
    path(query_fasta)
    path(search_db)
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