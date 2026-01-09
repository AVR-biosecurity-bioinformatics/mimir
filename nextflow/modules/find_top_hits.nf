process FIND_TOP_HITS {
    def module_name = "find_top_hits"
    // tag "-"
    // container "nanozoo/mmseqs2:14.7e284--11077ba"
    container "staphb/mmseqs2:18-8cc5c"

    input:
    path(query_fasta)
    path(target_fasta)
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
        ${task.memory.getMega()} \
        "${query_fasta}" \
        "${target_fasta}" \
        ${n_top_hits}
        
    """

}