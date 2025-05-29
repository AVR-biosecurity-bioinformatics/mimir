/*
Import and format internal (user-supplied) sequences
*/


//// modules to import
include { IMPORT_INTERNAL                                            } from '../modules/import_internal'


workflow GET_INTERNAL {

    take:

    ch_internal_seqs


    main:

    /*
    Import and format internal sequences
    */

    //// import internal sequences from file and check format
    IMPORT_INTERNAL (
        ch_internal_seqs
    )
    
    //// populate and chunk internal channel
    IMPORT_INTERNAL.out.fasta
        .splitText( by: params.input_chunk_size, file: true )
        .set { ch_internal_fasta }

    //// get internal sequence names for preferencing
    ch_internal_fasta
        .splitFasta( record: [ header:true ] )
        .map { record -> record.header }
        .collectFile ( name: 'internal_names.txt', newLine: true )
        .set { ch_internal_names }

    //// fill internal names channel if no internal sequences exist
    ch_internal_names = ch_internal_names.ifEmpty('no_file').first()
    
    
    emit:

    ch_internal_fasta
    ch_internal_names

}