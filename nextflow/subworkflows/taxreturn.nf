/*
 *  Run classic R::taxreturn pipeline
 */


//// modules to import
include { GET_NCBI_TAXONOMY                                         } from '../modules/get_ncbi_taxonomy'
// include { CHUNK_TAXON                                               } from '../modules/error_model'
include { FETCH_GENBANK                                             } from '../modules/fetch_genbank'
// include { FETCH_BOLD                                                } from '../modules/error_model'
// include { FETCH_MITO                                                } from '../modules/error_model'
include { TRIM_PHMM                                                 } from '../modules/trim_phmm'
include { FILTER_PHMM                                                 } from '../modules/filter_phmm'
include { FILTER_STOP                                                 } from '../modules/filter_stop'
include { COMBINE_CHUNKS                                                 } from '../modules/combine_chunks'
include { RESOLVE_SYNONYMS                                                 } from '../modules/resolve_synonyms'
include { REMOVE_CONTAM                                                 } from '../modules/remove_contam'



workflow TAXRETURN {

    take:

    seqs_internal


    main:

    //// dummy channel of taxon to test
    ch_taxon = Channel.from(
        "Scaptodrosophila",
        "Phylloxeridae"
        )

    //// make empty channels
    ch_genbank          = Channel.empty()
    ch_bold             = Channel.empty()
    ch_mito             = Channel.empty()
    ch_internal         = Channel.empty()

    //// parse file path parameters as channels
    ch_phmm = channel.fromPath( params.phmm_model, checkIfExists: true).first()

    //// get NCBI taxonomy file
    GET_NCBI_TAXONOMY ()

    /// use taxize to get taxonomic ranks from NCBI and BOLD: https://docs.ropensci.org/taxize/

    // //// split input taxon into units for parallelisation
    // CHUNK_TAXON (

    // )

    //// fetch genbank sequences for each chunk 
    FETCH_GENBANK (
        ch_taxon,
        GET_NCBI_TAXONOMY.out.db_file
    )

    ch_genbank = 
        ch_genbank 
        .concat ( FETCH_GENBANK.out.seqs )

    // //// fetch BOLD sequences for each chunk
    // FETCH_BOLD (

    // )

    // //// fetch mitochondrial genomes from NCBI for each chunk
    // FETCH_MITO (

    // )

    //// add internal sequences


    //// concat output channels into a single channel for PHMM alignment
    ch_genbank
        .concat ( ch_bold )
        .concat ( ch_mito )
        .concat ( ch_internal )
        .set { ch_raw_seqs }        


    //// trim model to primer sequences
    if ( params.trim_phmm ) {
        //// throw error if either primer sequence is not supplied
        if ( !params.primer_fwd || !params.primer_rev ) {
            println "*** params.primer_fwd = '$params.primer_fwd'; params.primer_rev = '$params.primer_rev' ***"
            error "*** ERROR: Both primer sequences must be given if '--trim_phmm' is set to 'true' ***"
        }
        
        TRIM_PHMM (
            ch_phmm, 
            params.primer_fwd,
            params.primer_rev,
            params.remove_primers
        )        
    }

    //// filter sequences in each chunk using PHMM model
    FILTER_PHMM (
        ch_raw_seqs,
        ch_phmm
    )

    //// optional: filter for stop codons
    if ( params.coding ){

        FILTER_STOP (
            FILTER_PHMM.out.seqs,
            "SGC4"
        )

        //// create chunks channel from FILTER_STOP output
        FILTER_STOP.out.seqs
            .map { taxon, index, type, seqs_file -> seqs_file }
            .collect()
            .set { ch_chunks }

    } else {
        
        //// create chunks channel from FILTER_PHMM output
        FILTER_PHMM.out.seqs
            .map { taxon, index, type, seqs_file -> seqs_file }
            .collect()
            .set { ch_chunks }
    }

    //// combine chunks together
    COMBINE_CHUNKS ( 
        ch_chunks 
    )
    
    //// resolve taxonomic synonyms
    RESOLVE_SYNONYMS ( 
        COMBINE_CHUNKS.out.seqs,
        GET_NCBI_TAXONOMY.out.db_path
    )

    //// remove contaminating sequences
    REMOVE_CONTAM (
        RESOLVE_SYNONYMS.out.seqs,
        GET_NCBI_TAXONOMY.out.db_file
    )

    //// prune large groups
    // PRUNE_GROUPS

    //// reformat names using taxonomic hierarchy
    // REFORMAT_NAMES

    //// summarise number of taxa in database
    // TAXA_SUMMARY

    //// train IDTAXA model
    // TRAIN_IDTAXA

    emit:

    ncbi_taxonomy = GET_NCBI_TAXONOMY.out.db_file


}