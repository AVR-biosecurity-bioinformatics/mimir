/*
Get sequences from BOLD
*/


//// modules to import
include { EXTRACT_BOLD                                               } from '../modules/extract_bold'
include { GET_BOLD_DATABASE                                          } from '../modules/get_bold_database'
include { MATCH_BOLD                                                 } from '../modules/match_bold'
include { MERGE_BOLD                                                 } from '../modules/merge_bold'


workflow GET_BOLD {

    take:

    ch_bold_db_path
    ch_bold_db_url
    ch_bold_targets
    ch_bold_query
    ch_lineageparents
    ch_synonyms


    main:

    /*
    Getting BOLD sequences and matching them to the NCBI database
    */
        
    //// get BOLD database files
    GET_BOLD_DATABASE (
        ch_bold_db_path,
        ch_bold_db_url
    )

    //// split BOLD database .tsv into chunks for processing
    GET_BOLD_DATABASE.out.bold_db
        .splitText (
            by: 500000, // TODO: make this a pipeline parameter
            keepHeader: true, 
            elem: 0,
            file: true
        )
        .set { ch_bold_db_chunks }

    //// combine each chunk with the taxon name and rank tibble from PARSE_TARGETS
    ch_bold_db_chunks
        .combine ( ch_bold_targets ) 
        .set { ch_extract_bold_input }

    //// extract sequences from BOLD database file
    EXTRACT_BOLD (
        ch_extract_bold_input,
        ch_bold_query,
        params.bold_idmethod_filter,
        params.min_length_input,
        params.max_length_input
    )

    //// match BOLD taxon names to NCBI taxon names
    MATCH_BOLD (
        EXTRACT_BOLD.out.tibble, 
        ch_lineageparents,
        ch_synonyms,
        params.placeholder_as_unclassified,
        params.digits_as_unclassified
    )

    //// refactor output channels so caching is valid 
    MATCH_BOLD.out.matching_data
        .flatten()
        .branch{ file ->
            fasta: file.fileName.name =~ /^bold_seqs/
            matching_taxids: file.fileName.name =~ /^matching_taxids/
            synchanges: file.fileName.name =~ /^synchanges/
        }
        .set{ ch_matched_bold }

    //// merge BOLD chunks into single .fasta and .csv files
    MERGE_BOLD (
        ch_matched_bold.fasta.collect(sort: true),
        ch_matched_bold.matching_taxids.collect(sort: true),
        ch_matched_bold.synchanges.collect(sort: true)
    )

    //// chunk BOLD sequences into smaller .fasta files for processing
    MERGE_BOLD.out.fasta
        .splitFasta( by: params.input_chunk_size, file: true )
        .set { ch_bold_fasta }

    
    emit:

    ch_bold_fasta
        

}