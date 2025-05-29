/*
Summarise database
*/


//// modules to import
include { CHECK_KEY_SPECIES                                          } from '../modules/check_key_species'
include { FASTA_TO_TABLE as CONCAT_FATES                             } from '../modules/fasta_to_table'
include { FASTA_TO_TABLE as CONCAT_SOURCES                           } from '../modules/fasta_to_table'
include { JOIN_SOURCES_FATES                                         } from '../modules/join_sources_fates'
include { SEQUENCE_TRACKER                                           } from '../modules/sequence_tracker'
include { VALIDATE_KEY_SPECIES                                       } from '../modules/validate_key_species'


workflow SUMMARISE_DATABASE {

    take:

    ch_final_database
    ch_sources
    ch_fates
    ch_key_species_list


    main:

    //// convert source .fasta files into a single .csv
    CONCAT_SOURCES (
        ch_sources,
        "sources"
    )

    //// concat sequence fates into one .csv
    CONCAT_FATES (
        ch_fates,
        "fates"
    )

    //// join fates and sources together into a single .csv
    JOIN_SOURCES_FATES (
        CONCAT_SOURCES.out.tsv,
        CONCAT_FATES.out.tsv
    )

    //// process the source and fate of every sequences through the pipeline
    SEQUENCE_TRACKER (
        JOIN_SOURCES_FATES.out.tsv
    )
  
    if ( params.key_species_list ){

        //// check for sequences belonging to key species
        CHECK_KEY_SPECIES (
            SEQUENCE_TRACKER.out.sf_meta,
            ch_key_species_list
        )

        //// validate sequences from key species
        VALIDATE_KEY_SPECIES (
            CHECK_KEY_SPECIES.out.seq_names.flatten(),
            ch_final_database,
            params.add_root
        )
        
    }


    emit:

    ch_out = SEQUENCE_TRACKER.out.taxa_summary

}