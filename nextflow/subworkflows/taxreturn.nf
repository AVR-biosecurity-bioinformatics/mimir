/*
 *  Run classic R::taxreturn pipeline
 */


//// modules to import
include { CHUNK_TAXON                                               } from '../modules/chunk_taxon'
include { COMBINE_CHUNKS                                                 } from '../modules/combine_chunks'
include { COUNT_GENBANK                                                 } from '../modules/count_genbank'
include { EXTRACT_BOLD                                                 } from '../modules/extract_bold'
// include { FETCH_BOLD                                                } from '../modules/fetch_bold'
include { FETCH_GENBANK                                             } from '../modules/fetch_genbank'
// include { FETCH_MITO                                                } from '../modules/fetch_mito'
include { FILTER_PHMM                                                 } from '../modules/filter_phmm'
include { FILTER_STOP                                                 } from '../modules/filter_stop'
include { GET_BOLD_DATABASE                                         } from '../modules/get_bold_database'
include { GET_NCBI_TAXONOMY                                         } from '../modules/get_ncbi_taxonomy'
include { MATCH_BOLD                                                 } from '../modules/match_bold'
include { MERGE_BOLD                                                 } from '../modules/merge_bold'
// include { PARSE_MARKER                                                 } from '../modules/parse_marker'
include { PARSE_TARGETS                                                 } from '../modules/parse_targets'
include { PRUNE_GROUPS                                                 } from '../modules/prune_groups'
include { REFORMAT_NAMES                                                 } from '../modules/reformat_names'
include { REMOVE_CONTAM                                                 } from '../modules/remove_contam'
include { REMOVE_EXACT_DUPLICATES                                                 } from '../modules/remove_exact_duplicates'
include { RESOLVE_SYNONYMS                                                 } from '../modules/resolve_synonyms'
include { TAXA_SUMMARY                                                 } from '../modules/taxa_summary'
include { TRAIN_IDTAXA                                                 } from '../modules/train_idtaxa'
include { TRIM_PHMM                                                 } from '../modules/trim_phmm'



workflow TAXRETURN {

    take:

    seqs_internal


    main:

    //// channel of taxon to test
    if ( params.target_taxa ) {
        // split comma_separated list into Groovy list
        ch_targets = params.target_taxa //// TODO: split comma-sep into list!
    } else {
        ch_targets = Channel.from(
            // "Drosophila", "Aphis"
            // "Scaptodrosophila", "Phylloxeridae"
            "Drosophilidae" 
        )
    }

    ch_target_ranks = params.target_ranks ?: Channel.of("no_ranks")

    //// ENTREZ key channel parsing
    ch_entrez_key = params.entrez_key ?: "no_key"

    //// BOLD database channel parsing
    if ( params.bold_db_path ) {
        ch_bold_db_path = channel.fromPath ( params.bold_db_path , checkIfExists: true, type: 'any' ).first()
    } else {
        ch_bold_db_path = "no_path"
    }

    ch_bold_db_url = params.bold_db_url ?: "no_url"

    //// make empty channels
    ch_genbank          = Channel.empty()
    ch_bold             = Channel.empty()
    ch_mito             = Channel.empty()
    ch_internal         = Channel.empty()

    //// parse file path parameters as channels
    ch_phmm             = channel.fromPath( params.phmm_model, checkIfExists: true ).first()

    //// get NCBI taxonomy file
    GET_NCBI_TAXONOMY ()

    //// convert input taxon/taxa into NCBI and BOLD tax IDs
    PARSE_TARGETS (
        ch_targets,
        ch_target_ranks,
        ch_entrez_key,
        GET_NCBI_TAXONOMY.out.synonyms
    )

    //// TODO: How does current PARSE_TARGETS implementation work with --target_taxa as a list?

    //// convert file text to channel value
    PARSE_TARGETS.out.taxon_name
        .map { name, rank -> 
            [ name.text, rank ] }
        .set { ch_taxon_namerank }



    // //// parse marker gene into formats understandable for each database
    // PARSE_MARKER (

    // )

    // //// query GenBank to get list of nucleotide IDs
    // COUNT_ENTREZ (

    // )


    // //// split input taxon into units for parallelisation
    // CHUNK_TAXON (
    //     ch_targets, 
    //     ch_entrez_key,
    //     GET_NCBI_TAXONOMY.out.db_path
    // )

    // //// split taxon list into a channel
    // CHUNK_TAXON.out.tax_list
    //     .splitText( by: 1 )
    //     .map { string -> string.trim() } // remove newline characters
    //     .set { ch_tax_chunks }

    if ( params.bold_db_path || params.bold_db_url ) {
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
            // .first() /// for testing only
            .set { ch_bold_db_chunks }

        //// extract sequences from BOLD database file
        EXTRACT_BOLD (
            ch_bold_db_chunks,
            PARSE_TARGETS.out.bold_names,
            PARSE_TARGETS.out.bold_rank
            // ch_marker
        )

        //// match BOLD taxon names to NCBI taxon names
        MATCH_BOLD (
            EXTRACT_BOLD.out.tibble, 
            GET_NCBI_TAXONOMY.out.lineageparents,
            GET_NCBI_TAXONOMY.out.synonyms
        )

        //// merge BOLD chunks into single .fasta and .csv files
        MERGE_BOLD (
            MATCH_BOLD.out.fasta.collect(),
            MATCH_BOLD.out.matching_taxids.collect(),
            MATCH_BOLD.out.synchanges.collect()
        )

        // ch_bold_seqs = EXTRACT_BOLD.out.tibble

    } else {
        //// extract sequences by querying BOLD API
        // FETCH_BOLD (

        // )
        //// TODO: Need to make sure this outputs same tibble format as EXTRACT_BOLD
    }

    //// combine BOLD sequences or dataframes into a single object

    


    // //// fetch count of genbank sequences for each chunk
    // COUNT_GENBANK (
    //     ch_tax_chunks
    // )

    // //// filter out chunks that return no sequences
    // COUNT_GENBANK.out.chunks_counts
    //     .map { taxon, count_file -> 
    //         int seq_count = count_file.getBaseName() as int
    //         [ taxon, count_file, seq_count ] }
    //     .filter { it[2] > 0 } // remove chunks with a sequence count of 0
    //     .map { taxon, count_file, seq_count -> [ taxon, count_file ] }
    //     .set { ch_filtered_chunks }

    // //// fetch genbank sequences for each chunk 
    // FETCH_GENBANK (
    //     ch_filtered_chunks,
    //     GET_NCBI_TAXONOMY.out.rankedlineage
    // )

    // //// split .fasta into subchunks for processing
    // ch_genbank 
    //     .concat ( FETCH_GENBANK.out.fasta )
    //     .splitFasta ( 
    //         by: params.subchunk_size,
    //         elem: 3,
    //         file: true
    //     )
    //     .set { ch_genbank }

    // // //// fetch BOLD sequences for each chunk
    // // FETCH_BOLD (

    // // )

    // // //// fetch mitochondrial genomes from NCBI for each chunk
    // // FETCH_MITO (

    // // )

    // //// add internal sequences


    // //// concat output channels into a single channel for PHMM alignment
    // ch_genbank
    //     .concat ( ch_bold )
    //     .concat ( ch_mito )
    //     .concat ( ch_internal )
    //     .set { ch_raw_seqs }        


    // //// trim model to primer sequences
    // if ( params.trim_phmm ) {
    //     //// throw error if either primer sequence is not supplied
    //     if ( !params.primer_fwd || !params.primer_rev ) {
    //         println "*** params.primer_fwd = '$params.primer_fwd'; params.primer_rev = '$params.primer_rev' ***"
    //         error "*** ERROR: Both primer sequences must be given if '--trim_phmm' is set to 'true' ***"
    //     }
        
    //     TRIM_PHMM (
    //         ch_phmm, 
    //         params.primer_fwd,
    //         params.primer_rev,
    //         params.remove_primers
    //     )        
    // }

    // //// filter sequences in each chunk using PHMM model
    // FILTER_PHMM (
    //     ch_raw_seqs,
    //     ch_phmm
    // )

    // //// optional: filter for stop codons
    // if ( params.coding && params.genetic_code ){

    //     FILTER_STOP (
    //         FILTER_PHMM.out.seqs
    //     )

    //     ch_filter_output = FILTER_STOP.out.seqs

    // } else {

    //     ch_filter_output = FILTER_PHMM.out.seqs
    
    // }

    // //// resolve taxonomic synonyms
    // RESOLVE_SYNONYMS ( 
    //     FILTER_PHMM.out.seqs,
    //     GET_NCBI_TAXONOMY.out.db_path
    // )

    // //// collect chunks into a list
    // RESOLVE_SYNONYMS.out.seqs
    //     .collect()
    //     .set { ch_chunks }

    // //// combine sequence chunks together
    // COMBINE_CHUNKS ( 
    //     ch_chunks
    // )
 
    // //// remove exact duplicate sequences if they exist
    // REMOVE_EXACT_DUPLICATES (
    //     COMBINE_CHUNKS.out.fasta
    // )

    // //// remove contaminating sequences
    // REMOVE_CONTAM (
    //     REMOVE_EXACT_DUPLICATES.out.fasta,
    //     GET_NCBI_TAXONOMY.out.rankedlineage
    // )

    // //// prune large groups
    // PRUNE_GROUPS (
    //     REMOVE_CONTAM.out.seqs
    // )

    // //// reformat names using taxonomic hierarchy
    // REFORMAT_NAMES (
    //     PRUNE_GROUPS.out.seqs,
    //     GET_NCBI_TAXONOMY.out.rankedlineage
    // )

    // //// summarise number of taxa in database
    // TAXA_SUMMARY (
    //     REFORMAT_NAMES.out.seqs
    // )

    // //// train IDTAXA model
    // if ( params.train_idtaxa ) {
    //     TRAIN_IDTAXA (
    //         REFORMAT_NAMES.out.seqs,
    //         GET_NCBI_TAXONOMY.out.rankedlineage
    //     )
    // }
    
    // //// TODO: add proper model channel handling (with conditions)
    // ch_models = 
    //     Channel.empty()
    

    emit:

    // bold_db = GET_BOLD_DATABASE.out.bold_db
    ncbi_taxonomy = GET_NCBI_TAXONOMY.out.rankedlineage
    // curated_fasta = REFORMAT_NAMES.out.fasta
    // taxa_summary = TAXA_SUMMARY.out.csv
    // idtaxa_model = TRAIN_IDTAXA.out.model

}