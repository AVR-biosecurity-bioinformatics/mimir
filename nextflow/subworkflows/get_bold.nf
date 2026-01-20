/*
Get sequences from BOLD
*/


//// modules to import
include { EXTRACT_BOLD                                               } from '../modules/extract_bold'
include { EXTRACT_BOLD_AWK                                           } from '../modules/extract_bold_awk'
include { GET_BOLD_DATABASE                                          } from '../modules/get_bold_database'
include { MATCH_BOLD                                                 } from '../modules/match_bold'


workflow GET_BOLD {

    take:

    ch_bold_db_path
    ch_bold_db_url
    ch_bold_targets
    ch_bold_query
    ch_lineageparents
    ch_filteredsynonyms
    ch_genbank_accessions

    main:

    /*
    Getting BOLD sequences and matching them to the NCBI database
    */
        
    //// get BOLD database files
    GET_BOLD_DATABASE (
        ch_bold_db_path,
        ch_bold_db_url
    )

    //// extract sequences from whole BOLD database file using awk
    EXTRACT_BOLD_AWK (
        GET_BOLD_DATABASE.out.bold_db,
        ch_bold_targets,
        ch_bold_query,
        params.min_length_input,
        params.max_length_input
    )

    //// split into 1m record chunks 
    EXTRACT_BOLD_AWK.out.tsv
        .filter { it.size() > 0 }
        .splitText( by: 1000000, keepHeader: true, file: true )
        .set { ch_bold_db_extracted }

    //// process BOLD sequences, remove duplicated GenBank sequences (if needed), and harmonise to NCBI taxonomy
    MATCH_BOLD (
        ch_bold_db_extracted, 
        ch_lineageparents.first(),
        ch_filteredsynonyms.first(),
        ch_genbank_accessions.first(),
        params.placeholder_as_unclassified,
        params.digits_as_unclassified,
        params.bold_idmethod_filter
    )

    MATCH_BOLD.out.fasta
        .collectFile ( name: 'bold_seqs.fasta' )
        .set { ch_bold_fasta }
    
    MATCH_BOLD.out.gb_dups
        .collectFile ( name: 'gb_dups.csv', keepHeader: true )
        .set { ch_gb_dups }

    //// merge other tables
    MATCH_BOLD.out.matching_taxids
        .collectFile ( 
            name: "bold_matches.csv",
            storeDir: "./output/results",
            keepHeader: true, 
            skip: 1
        )

    MATCH_BOLD.out.synchanges
        .collectFile ( 
            name: "bold_synchanges.csv",
            storeDir: "./output/results",
            keepHeader: true, 
            skip: 1
        )
   
    emit:

    ch_bold_fasta
    ch_gb_dups

}