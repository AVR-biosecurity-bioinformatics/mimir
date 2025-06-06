process MERGE_BOLD {
    def module_name = "merge_bold"
    // tag "-"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(fasta_list, name: 'seqs*.fasta')
    path(matching_taxids_list, name: 'matching_taxids*.csv')
    path(synchanges_list, name: 'synchanges*.csv')

    output: 
    path("bold_seqs.merged.fasta"),                  emit: fasta
    path("matching_taxids.merged.csv"),              emit: matching_taxids
    path("synchanges.merged.csv"),                   emit: synchanges

    // publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env bash
    
    ### run module code
    bash ${module_name}.sh \
        ${projectDir} \
        ${task.cpus} \
        "${fasta_list}" \
        "${matching_taxids_list}" \
        "${synchanges_list}"

    """
}