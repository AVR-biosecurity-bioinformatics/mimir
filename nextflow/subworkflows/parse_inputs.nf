/*
Parse pipeline inputs and fetch taxonomic and marker-related information
*/


//// modules to import
include { ALIGN_PRIMERS_TO_SEED                                      } from '../modules/align_primers_to_seed'
include { BUILD_TRIMMED_HMM                                          } from '../modules/build_trimmed_hmm'
include { GET_NCBI_TAXONOMY                                          } from '../modules/get_ncbi_taxonomy'
include { MAKE_GENCODES                                              } from '../modules/make_gencodes'
include { MAKE_LINEAGEPARENTS                                        } from '../modules/make_lineageparents'
include { MAKE_NAMES                                                 } from '../modules/make_names'
include { MAKE_NODES                                                 } from '../modules/make_nodes'
include { MAKE_RANKEDLINEAGE_NONAME                                  } from '../modules/make_rankedlineage_noname'
include { MAKE_RANKEDLINEAGE                                         } from '../modules/make_rankedlineage'
include { MAKE_SYNONYMS                                              } from '../modules/make_synonyms'
include { MAKE_TAXIDNAMERANK                                         } from '../modules/make_taxidnamerank'
include { MAKE_TAXIDNAMES                                            } from '../modules/make_taxidnames'
include { MATCH_BOLD                                                 } from '../modules/match_bold'
include { PARSE_MARKER                                               } from '../modules/parse_marker'
include { PARSE_TARGETS                                              } from '../modules/parse_targets'
include { PROCESS_PRIMERS                                            } from '../modules/process_primers'
include { STOCKHOLM_TO_FASTA                                         } from '../modules/stockholm_to_fasta'
include { TRIM_HMM_SEED                                              } from '../modules/trim_hmm_seed'


workflow PARSE_INPUTS {

    take:

    ch_dummy
    ch_entrez_key
    ch_targets
    ch_marker

    main:

    /*
    Process NCBI Taxonomy data
    */

    //// get NCBI taxonomy file
    GET_NCBI_TAXONOMY (
        ch_dummy
    )

    //// make rankedlineage tibble
    MAKE_RANKEDLINEAGE (
        GET_NCBI_TAXONOMY.out.db_path
    )

    //// make nodes tibble
    MAKE_NODES (
        GET_NCBI_TAXONOMY.out.db_path
    )

    //// make taxidnames tibble
    MAKE_TAXIDNAMES (
        GET_NCBI_TAXONOMY.out.db_path
    )

    //// make names tibble
    MAKE_NAMES (
        GET_NCBI_TAXONOMY.out.db_path
    )

    //// make gencodes tibble
    MAKE_GENCODES (
        MAKE_RANKEDLINEAGE.out.rankedlineage,
        MAKE_NODES.out.nodes
    )

    //// make taxidnamerank tibble
    MAKE_TAXIDNAMERANK (
        MAKE_NODES.out.nodes,
        MAKE_TAXIDNAMES.out.taxidnames
    )

    //// make synonyms tibble
    MAKE_SYNONYMS (
        MAKE_TAXIDNAMERANK.out.taxidnamerank,
        MAKE_NAMES.out.names
    )

    //// make lineageparents tibble
    MAKE_LINEAGEPARENTS (
        MAKE_RANKEDLINEAGE.out.rankedlineage,
        MAKE_TAXIDNAMERANK.out.taxidnamerank
    )

    //// make rankedlineage_noname tibble
    MAKE_RANKEDLINEAGE_NONAME (
        MAKE_RANKEDLINEAGE.out.rankedlineage,
        MAKE_TAXIDNAMERANK.out.taxidnamerank
    )

    /* 
    Parse inputs
    */

    //// convert input taxon/taxa into NCBI and BOLD tax IDs
    PARSE_TARGETS (
        ch_targets,
        ch_entrez_key,
        MAKE_SYNONYMS.out.synonyms,
        MAKE_NODES.out.nodes
    )

    //// convert file text to channel value
    PARSE_TARGETS.out.taxon_id_rank
        .map { id, rank -> 
            [ id.text.trim(), rank ] }
        .set { ch_taxon_idrank }

    //// parse marker gene into formats understandable for each database/process
    PARSE_MARKER (
        ch_marker
    )

    //// collect target gencodes .csvs into a single file
    PARSE_TARGETS.out.gencodes
        .collectFile ( name: 'gencodes.csv', keepHeader: true, skip: 1 )
        .set { ch_target_gencodes }

    //// processes relevant to trimming sequences
    if ( params.trim_to_primers ) {
        
        //// get disambiguated and translated primer sequences
        PROCESS_PRIMERS (
            params.primer_fwd,
            params.primer_rev,
            ch_target_gencodes,
            PARSE_MARKER.out.type
        )

        //// convert seed alignment to fasta format
        STOCKHOLM_TO_FASTA (
            PARSE_MARKER.out.seed
        )

        //// align translated primers to PHMM seed alignment
        ALIGN_PRIMERS_TO_SEED (
            STOCKHOLM_TO_FASTA.out.fasta,
            PROCESS_PRIMERS.out.translated
        )

        //// trim seed alignment to primer region
        TRIM_HMM_SEED (
            ALIGN_PRIMERS_TO_SEED.out.fasta
        )

        //// build trimmed PHMM from trimmed seed alignment
        BUILD_TRIMMED_HMM (
            TRIM_HMM_SEED.out.fasta
        )

        ch_trimmed_phmm = BUILD_TRIMMED_HMM.out.hmm.first()

    } else {
        ch_trimmed_phmm = Channel.empty()
    }

    emit:

    ch_taxon_idrank                 
    ch_bold_targets                 = PARSE_TARGETS.out.bold
    ch_genbank_query                = PARSE_MARKER.out.genbank_query
    ch_bold_query                   = PARSE_MARKER.out.bold_query
    ch_marker_coding                = PARSE_MARKER.out.coding
    ch_marker_type                  = PARSE_MARKER.out.type
    ch_rankedlineage_noname         = MAKE_RANKEDLINEAGE_NONAME.out.rankedlineage_noname
    ch_lineageparents               = MAKE_LINEAGEPARENTS.out.lineageparents
    ch_synonyms                     = MAKE_SYNONYMS.out.synonyms
    ch_gencodes                     = MAKE_GENCODES.out.gencodes
    ch_full_phmm                    = PARSE_MARKER.out.phmm
    ch_trimmed_phmm

}