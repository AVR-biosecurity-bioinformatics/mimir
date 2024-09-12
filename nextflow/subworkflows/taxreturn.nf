/*
 *  Run classic R::taxreturn pipeline
 */


//// modules to import
include { GET_NCBI_TAXONOMY                                         } from '../modules/get_ncbi_taxonomy'
// include { CHUNK_TAXON                                               } from '../modules/error_model'
include { FETCH_GENBANK                                             } from '../modules/fetch_genbank'
// include { FETCH_BOLD                                                } from '../modules/error_model'
// include { FETCH_MITO                                                } from '../modules/error_model'


workflow TAXRETURN {

    take:

    seqs_internal


    main:

    //// dummy channel of taxon to test
    ch_taxon = Channel.from("Scaptodrosophila","Carpophilus")


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

    // //// fetch BOLD sequences for each chunk
    // FETCH_BOLD (

    // )

    // //// fetch mitochondrial genomes from NCBI for each chunk
    // FETCH_MITO (

    // )

    //// add internal sequences


    //// concat output channels into a single channel for PHMM alignment


    //// optional: trim model to primer sequences


    //// align sequences in each chunk to PHMM model


    //// optional: filter for stop codons


    //// combine chunks together

    
    //// resolve synonyms


    //// remove contaminating sequences


    //// prune large groups


    //// reformat names using taxonomic hierarchy


    //// summarise number of taxa in database


    //// train IDTAXA model


    emit:

    ncbi_taxonomy = GET_NCBI_TAXONOMY.out.db_file


}