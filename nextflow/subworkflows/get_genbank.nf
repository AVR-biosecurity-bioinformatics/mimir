/*
Get sequences from GenBank 
*/


//// modules to import
include { FETCH_GENBANK                                              } from '../modules/fetch_genbank'
include { QUERY_GENBANK                                              } from '../modules/query_genbank'
include { RENAME_GENBANK                                             } from '../modules/rename_genbank'


workflow GET_GENBANK {

    take:

    ch_taxon_idrank
    ch_genbank_query
    ch_entrez_key
    ch_rankedlineage_noname


    main:

    /*
    Get Genbank sequences
    */

    //// query GenBank to get list of nucleotide IDs (including mitochondrial genomes if requested)
    QUERY_GENBANK (
        ch_taxon_idrank,
        ch_genbank_query,
        params.min_length_input,
        params.max_length_input,
        params.use_mito
    )

    //// split list of accessions for fetching in chunks
    QUERY_GENBANK.out.seq_acc
        .splitText( by: params.genbank_fetch_size, file: true )
        .set { ch_genbank_acc_chunks }

    //// fetch Genbank sequences as .fasta + taxid list
    FETCH_GENBANK (
        ch_genbank_acc_chunks,
        ch_entrez_key
    )

    //// reformat sequence names to contain taxonomic lineage
    RENAME_GENBANK (
        FETCH_GENBANK.out.fetched_seqs,
        ch_rankedlineage_noname,
        params.placeholder_as_unclassified,
        params.digits_as_unclassified
    )

    
    emit:

    ch_genbank_fasta                = RENAME_GENBANK.out.fasta
    

}