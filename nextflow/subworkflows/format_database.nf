/*
Format final database
*/


//// modules to import
include { ALIGN_DATABASE_CORE                                        } from '../modules/align_database_core'
include { ALIGN_DATABASE_OTHER                                       } from '../modules/align_database_other'
include { COMBINE_CHUNKS                                             } from '../modules/combine_chunks'
include { FORMAT_OUTPUT                                              } from '../modules/format_output'
include { GET_DATABASE_CORE                                          } from '../modules/get_database_core'
include { TRAIN_IDTAXA                                               } from '../modules/train_idtaxa'


workflow FORMAT_DATABASE {

    take:

    ch_final_seqs


    main:

    //// collect final sequence .fasta into a single element
    ch_final_seqs
        .flatten()
        .collect()
        .set { ch_combine_input }

    //// combine species-level .fasta files into a single file
    COMBINE_CHUNKS ( 
        ch_combine_input,
        "false"
    )

    /// align databse if requested
    if ( params.aligned_output ) {
        
        //// whole-database alignment method based on mafft-sparsecore.rb (adding sequences to a core alignment)

        //// define 'core' and 'other' sequences based on taxonomy and length
        GET_DATABASE_CORE (
            COMBINE_CHUNKS.out.fasta
        )

        //// align core sequences
        ALIGN_DATABASE_CORE (
            GET_DATABASE_CORE.out.core
        )

        //// add other sequences to core alignment
        ALIGN_DATABASE_OTHER (
            ALIGN_DATABASE_CORE.out.fasta,
            GET_CORE_SEQUENCES.out.other
        )

        ALIGN_DATABASE_OTHER.out.fasta
            .set { ch_aligned_database }

        //// save aligned database
        if ( params.save_intermediate ) {
            ch_aligned_database
                .collectFile ( 
                    name: "aligned_database.fasta",
                    storeDir: "./output/results"
                )
        }

        ch_formatting_input = ch_aligned_database

    } else {

        ch_formatting_input = COMBINE_CHUNKS.out.fasta
    }

    //// format database output
    FORMAT_OUTPUT (
        ch_formatting_input,
        params.add_root,
        params.aligned_output,
        params.compressed_output
    )

    //// save final output
    FORMAT_OUTPUT.out.fasta
        .collectFile ( 
            name: "format_output.fasta",
            storeDir: "./output/results"
        )

    //// train IDTAXA model
    if ( params.train_idtaxa ) {
    
        TRAIN_IDTAXA (
            FORMAT_OUTPUT.out.fasta,
            params.idtaxa_max_group_size,
            params.idtaxa_max_iterations,
            params.idtaxa_allow_group_removal
        )

        ch_idtaxa_model = TRAIN_IDTAXA.out.model
    
    } else {

        ch_idtaxa_model = Channel.empty()

    }
    

    emit:

    ch_final_database = FORMAT_OUTPUT.out.fasta
    ch_idtaxa_model

}