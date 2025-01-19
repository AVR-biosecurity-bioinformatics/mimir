/*
 *  Run classic R::taxreturn pipeline
 */


//// modules to import
include { ADD_TAXID_GENBANK                                                 } from '../modules/add_taxid_genbank'
include { ALIGN_BATCH as ALIGN_SPECIES                                                 } from '../modules/align_batch'
include { ALIGN_SINGLE as ALIGN_OUTPUT                                                 } from '../modules/align_single'
include { CLUSTER_SEQUENCES                                                 } from '../modules/cluster_sequences'
include { COMBINE_CHUNKS as COMBINE_CHUNKS_1                            } from '../modules/combine_chunks'
include { COMBINE_CHUNKS as COMBINE_CHUNKS_2                           } from '../modules/combine_chunks'
include { EXTRACT_BOLD                                                 } from '../modules/extract_bold'
// include { FETCH_BOLD                                                } from '../modules/fetch_bold'
include { FETCH_GENBANK                                             } from '../modules/fetch_genbank'
// include { FETCH_MITO                                                } from '../modules/fetch_mito'
include { FILTER_HMM                                                 } from '../modules/filter_hmm'
// include { FILTER_PHMM                                                 } from '../modules/filter_phmm'
// include { FILTER_STOP                                                 } from '../modules/filter_stop'
include { FORMAT_OUTPUT                                         } from '../modules/format_output'
include { GET_BOLD_DATABASE                                         } from '../modules/get_bold_database'
include { GET_NCBI_TAXONOMY                                         } from '../modules/get_ncbi_taxonomy'
include { HMMSEARCH                                         } from '../modules/hmmsearch'
include { IMPORT_INTERNAL                                                 } from '../modules/import_internal'
include { MATCH_BOLD                                                 } from '../modules/match_bold'
// include { MATCH_INTERNAL                                                 } from '../modules/match_internal'
include { MERGE_BOLD                                                 } from '../modules/merge_bold'
include { MERGE_SPLITS                                                 } from '../modules/merge_splits'
include { PARSE_MARKER                                                 } from '../modules/parse_marker'
include { PARSE_TARGETS                                                 } from '../modules/parse_targets'
include { PRUNE_GROUPS                                                 } from '../modules/prune_groups'
include { QUERY_GENBANK                                                 } from '../modules/query_genbank'
include { REFORMAT_NAMES                                                 } from '../modules/reformat_names'
include { REMOVE_AMBIGUOUS                                                 } from '../modules/remove_ambiguous'
include { REMOVE_EXACT_DUPLICATES                                                 } from '../modules/remove_exact_duplicates'
include { REMOVE_SEQ_OUTLIERS                                                 } from '../modules/remove_seq_outliers'
include { REMOVE_TAX_OUTLIERS                                                 } from '../modules/remove_tax_outliers'
include { REMOVE_UNCLASSIFIED                                                 } from '../modules/remove_unclassified'
include { RENAME_GENBANK                                                 } from '../modules/rename_genbank'
include { RESOLVE_SYNONYMS                                                 } from '../modules/resolve_synonyms'
include { SORT_BY_LINEAGE                                                } from '../modules/sort_by_lineage'
include { SPLIT_BY_RANK as SPLIT_BY_SPECIES                                                } from '../modules/split_by_rank'
include { SUMMARISE_COUNTS                                                 } from '../modules/summarise_counts'
include { SUMMARISE_TAXA                                                 } from '../modules/summarise_taxa'
include { TRAIN_IDTAXA                                                 } from '../modules/train_idtaxa'
include { TRANSLATE_SEQUENCES                                                 } from '../modules/translate_sequences'
// include { TRIM_PHMM                                                 } from '../modules/trim_phmm'



workflow TAXRETURN {

    take:

    seqs_internal


    main:

    /*
    Define channels
    */

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

    //// marker channel parsing
    if ( params.marker ) {
        ch_marker = params.marker
    } else {
        error "*** '--marker' must be specified ***"
    }

    //// ENTREZ key channel parsing
    ch_entrez_key = params.entrez_key ?: "no_key"

    //// BOLD database channel parsing
    if ( params.bold_db_path ) {
        ch_bold_db_path = channel.fromPath ( params.bold_db_path , checkIfExists: true, type: 'any' ).first()
    } else {
        ch_bold_db_path = "no_path"
    }

    ch_bold_db_url = params.bold_db_url ?: "no_url"

    //// internal sequence channel parsing
    if ( params.internal_seqs ) {
        ch_internal_seqs = channel.fromPath ( params.internal_seqs, checkIfExists: true, type: 'file').first()
    } else {
        ch_internal_seqs = Channel.empty()
    }

    //// make empty channels
    ch_genbank_fasta          = Channel.empty()
    ch_bold_fasta             = Channel.empty()
    ch_mito_fasta             = Channel.empty()
    ch_genome_fasta           = Channel.empty()
    ch_internal_fasta         = Channel.empty()
    ch_internal_names         = Channel.empty()
    ch_input_seqs             = Channel.empty()

    //// parse file path parameters as channels
    ch_phmm             = channel.fromPath( params.phmm_model, checkIfExists: true ).first()

    /*
    Set up pipeline
    */

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
    PARSE_TARGETS.out.taxon_id_rank
        .map { id, rank -> 
            [ id.text.trim(), rank ] }
        .set { ch_taxon_idrank }

    //// parse marker gene into formats understandable for each database/process
    PARSE_MARKER (
        ch_marker
    )

    /*
    Get Genbank sequences
    */
    if ( params.use_genbank ) {
        //// query GenBank to get list of nucleotide IDs (including mitochondrial genomes if requested)
        QUERY_GENBANK (
            ch_taxon_idrank,
            PARSE_MARKER.out.genbank_query
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

        //// add NCBI taxid to header of Genbank sequences 
        ADD_TAXID_GENBANK (
            FETCH_GENBANK.out.fetched_seqs
        )

        //// reformat sequence names to contain taxonomic lineage
        RENAME_GENBANK (
            ADD_TAXID_GENBANK.out.fasta,
            GET_NCBI_TAXONOMY.out.rankedlineage_noname
        )

        //// populate empty ch_genbank_fasta channel
        RENAME_GENBANK.out.fasta
            .set { ch_genbank_fasta }
    }
    

    //// count number of sequences downloaded from Genbank
    ch_count_genbank = ch_genbank_fasta.countFasta().combine(["genbank"])
    
    /*
    Getting BOLD sequences and matching them to the NCBI database
    */

    if ( ( params.bold_db_path || params.bold_db_url ) && params.use_bold ) {
        
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

        //// extract sequences from BOLD database file
        EXTRACT_BOLD (
            ch_bold_db_chunks,
            PARSE_TARGETS.out.bold_names,
            PARSE_TARGETS.out.bold_rank,
            PARSE_MARKER.out.bold_query
        )

        //// match BOLD taxon names to NCBI taxon names
        MATCH_BOLD (
            EXTRACT_BOLD.out.tibble, 
            GET_NCBI_TAXONOMY.out.lineageparents,
            GET_NCBI_TAXONOMY.out.synonyms
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

        

    } else {
        //// extract sequences by querying BOLD API
        // FETCH_BOLD (



        // )
        
    }

    //// count number of BOLD sequences
    ch_count_bold = ch_bold_fasta.countFasta().combine(["bold"])

    //// count number of mitochondrial sequences 
    ch_count_mito = ch_mito_fasta.countFasta().combine(["mito"])

    //// count number of genome assembly-derived sequences 
    ch_count_genome = ch_genome_fasta.countFasta().combine(["genome"])



    /*
    Import and format internal sequences
    */

    if ( ch_internal_seqs ){
        //// import internal sequences from file and check format
        IMPORT_INTERNAL (
            ch_internal_seqs
        )

    //     //// match the taxonomy of internal sequences to NCBI taxonomy
    //     MATCH_INTERNAL (
    //         IMPORT_INTERNAL.out.fasta
    //     )
        
        //// populate and chunk internal channel
        IMPORT_INTERNAL.out.fasta
            .splitText( by: params.input_chunk_size, file: true )
            .set { ch_internal_fasta }
    }
    
    //// get internal sequence names for preferencing
    ch_internal_fasta
        .splitFasta( record: [ header:true ] )
        .map { record -> record.header }
        .collectFile ( name: 'internal_names.txt', newLine: true )
        .set { ch_internal_names }

    //// fill internal names channel if no internal sequences exist
    ch_internal_names = ch_internal_names.ifEmpty('no_file').first()
    
    //// count number of internal sequences 
    ch_count_internal = ch_internal_fasta.countFasta().combine(["internal"])



    //// count number of external input sequences
    ch_count_external = ch_genbank_fasta
        .mix ( ch_bold_fasta )
        .mix ( ch_mito_fasta )
        .mix ( ch_genome_fasta )
        .countFasta()
        .combine(["external"])

    //// concat output channels into a single channel for PHMM alignment
    ch_input_seqs
        .mix ( ch_genbank_fasta )
        .mix ( ch_bold_fasta )
        .mix ( ch_mito_fasta )
        .mix ( ch_genome_fasta )
        .mix ( ch_internal_fasta )
        .set { ch_input_seqs }  

    //// count number of total input sequences
    ch_count_input = ch_input_seqs.countFasta().combine(["input"])      

    //// combine and save all unfiltered input sequences a single .fasta 
    if ( params.save_input || params.save_intermediate ) {
        ch_input_seqs
            .collectFile ( 
                name: "input_sequences.fasta",
                storeDir: "./output/results"
            )
    }

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

    ////
    ch_input_seqs
        .set { ch_remove_unclassified_input }

    //// remove unclassified sequences
    if ( params.remove_unclassified == "all_ranks" || params.remove_unclassified == "any_ranks" || params.remove_unclassified == "terminal" ) {
        REMOVE_UNCLASSIFIED (
            ch_remove_unclassified_input,
            params.remove_unclassified
        )        

        ch_count_remove_unclassified = REMOVE_UNCLASSIFIED.out.fasta.countFasta().combine(["remove_unclassified"])

        //// combine and save intermediate file 
        if ( params.save_intermediate ) {
            REMOVE_UNCLASSIFIED.out.fasta
                .collectFile ( 
                    name: "remove_unclassified.fasta",
                    storeDir: "./output/results",
                    cache: 'lenient'
                )
        }

        REMOVE_UNCLASSIFIED.out.fasta
            .set { ch_remove_ambiguous_input }

    } else {

        ch_count_remove_unclassified = Channel.of(["NA", "remove_unclassified"])

        ch_remove_unclassified_input
            .set { ch_remove_ambiguous_input }
    }
    //// remove ambiguous 
    if ( params.remove_ambiguous ){
        REMOVE_AMBIGUOUS (
            ch_remove_ambiguous_input
        )

        ch_count_remove_ambiguous = REMOVE_AMBIGUOUS.out.fasta.countFasta().combine(["remove_ambiguous"])

        //// combine and save intermediate file 
        if ( params.save_intermediate ) {
            REMOVE_AMBIGUOUS.out.fasta
                .collectFile ( 
                    name: "remove_ambiguous.fasta",
                    storeDir: "./output/results",
                    cache: 'lenient'
                )
        }

        REMOVE_AMBIGUOUS.out.fasta
            .set { ch_translate_sequences_input }

    } else {

        ch_count_remove_ambiguous = Channel.of(["NA", "remove_ambiguous"])

        ch_remove_ambiguous_input
            .set { ch_translate_sequences_input }
    }

    //// translate sequences in all six frames
    TRANSLATE_SEQUENCES (
        ch_translate_sequences_input,
        PARSE_MARKER.out.coding,
        PARSE_MARKER.out.type,
        GET_NCBI_TAXONOMY.out.ncbi_gencodes
    )

    //// search translated sequences for hits against PHMM model
    HMMSEARCH (
        TRANSLATE_SEQUENCES.out.translations,
        PARSE_MARKER.out.phmm
    )

    //// filter sequences by HMMSEARCH output
    FILTER_HMM (
        HMMSEARCH.out.hmmer_output,
        params.hmm_max_evalue,
        params.hmm_min_score,
        params.hmm_max_hits,
        params.hmm_min_acc,
        params.hmm_max_gap
    )

    //// combine and save intermediate file 
    if ( params.save_intermediate ) {
        FILTER_HMM.out.fasta
            .collectFile ( 
                name: "filter_hmm.fasta",
                storeDir: "./output/results",
                cache: 'lenient'
            )
    }

    //// combine and save excluded sequences
    if ( params.save_intermediate ) {
        FILTER_HMM.out.removed_fasta
            .collectFile ( 
                name: "filter_hmm.removed.fasta",
                storeDir: "./output/results",
                cache: 'lenient'
            )
    }

    //// combine and save excluded sequence table
    if ( params.save_intermediate ) {
        FILTER_HMM.out.removed_csv
            .collectFile ( 
                name: "filter_hmm.removed.csv",
                storeDir: "./output/results",
                cache: 'lenient',
                keepHeader: true, 
                skip: 1
            )
    }

    ////

    //// count number of sequences passing PHMM filter
    ch_count_filter_hmm = FILTER_HMM.out.fasta.countFasta().combine(["filter_hmm"]) 

    ch_filter_output = FILTER_HMM.out.fasta
        .filter{ it.size()>0 } // remove empty files
        .collect()

    //// combine sequences into one .fasta file and dealign
    COMBINE_CHUNKS_1 ( 
        ch_filter_output,
        "true"
    )
 
    //// remove exact duplicate sequences if they exist
    REMOVE_EXACT_DUPLICATES (
        COMBINE_CHUNKS_1.out.fasta
    )

    //// combine and save intermediate file 
    if ( params.save_intermediate ) {
        REMOVE_EXACT_DUPLICATES.out.fasta
            .collectFile ( 
                name: "remove_exact_duplicates.fasta",
                storeDir: "./output/results"
            )
    }

    //// count number of sequences passing exact deduplication
    ch_count_remove_exact = REMOVE_EXACT_DUPLICATES.out.fasta.countFasta().combine(["remove_exact"])

    //// cluster sequences into OTUs with mmseqs2
    CLUSTER_SEQUENCES (
        REMOVE_EXACT_DUPLICATES.out.fasta
    )

    //// remove taxonomic outliers from sequence clusters
    REMOVE_TAX_OUTLIERS (
        REMOVE_EXACT_DUPLICATES.out.fasta,
        CLUSTER_SEQUENCES.out.tsv,
        GET_NCBI_TAXONOMY.out.rankedlineage // NOTE: remove this channel and update process, as not needed
    )

    //// combine and save intermediate file 
    if ( params.save_intermediate ) {
        REMOVE_TAX_OUTLIERS.out.fasta
            .collectFile ( 
                name: "remove_tax_outliers.fasta",
                storeDir: "./output/results"
            )
    }

    //// count number of sequences passing taxonomic decontamination
    ch_count_remove_tax_outliers = REMOVE_TAX_OUTLIERS.out.fasta.countFasta().combine(["remove_tax_outliers"])

    // ch_split_species_input.take( 3 ).view()

    //// sort .fasta by lineage string to reduce the number of files that need to be merged after splitting 
    SORT_BY_LINEAGE (
        REMOVE_TAX_OUTLIERS.out.fasta
    )

    //// chunk input into SPLIT_BY_SPECIES to parallelise
    SORT_BY_LINEAGE.out.fasta
        .splitFasta ( by: params.split_rank_chunk, file: true )
        .set { ch_split_species_input }

    //// split .fasta by taxonomic lineage down to species level
    SPLIT_BY_SPECIES (
        ch_split_species_input,
        "species"
    )

    //// group files with the same name
    SPLIT_BY_SPECIES.out.fasta
        .flatten()
        .map { file ->
            [ file.name , file ] }
        .groupTuple( by: 0 )
        // branch channel depending on the number of files in each tuple
        .branch { file_name, file_list ->
            no_merge: file_list.size() == 1
                return file_list
            merge: file_list.size() > 1
                return file_list
        } 
        .set { ch_split_species_output }

    //// group merging files into groups of 1000 species for merging
    ch_split_species_output.merge
        .buffer ( size: 1000, remainder: true )
        .flatten()
        .collect()
        .set { ch_merge_splits_input }

    //// merge files from the same species across the different chunks
    MERGE_SPLITS (
        ch_merge_splits_input
    )

    //// group species-level .fasta files into batches for aligning and pruning
    ch_split_species_output.no_merge
        .flatten()
        .mix( MERGE_SPLITS.out.fasta.flatten() )
        .buffer( size: 100, remainder: true ) 
        .set { ch_align_input }

    //// align species-level .fasta in batches
    ALIGN_SPECIES (
        ch_align_input
    )

    //// remove sequence outliers from species clusters
    REMOVE_SEQ_OUTLIERS (
        ALIGN_SPECIES.out.aligned_fasta,
        params.dist_threshold
    )

    //// combine and save intermediate file 
    if ( params.save_intermediate ) {
        REMOVE_SEQ_OUTLIERS.out.retained_fasta
            .flatten()
            .collectFile ( 
                name: "remove_seq_outliers.fasta",
                storeDir: "./output/results"
            )
    }

    //// count number of sequences passing taxonomic decontamination
    ch_count_remove_seq_outliers = REMOVE_SEQ_OUTLIERS.out.retained_fasta.flatten().countFasta().combine(["remove_seq_outliers"])

    //// prune large groups, preferentially retaining internal sequences
    PRUNE_GROUPS (
        REMOVE_SEQ_OUTLIERS.out.retained_fasta,
        // ALIGN_SPECIES.out.aligned_fasta,
        ch_internal_names
    )

    //// save PRUNE_GROUPS output
    PRUNE_GROUPS.out.fasta
        .flatten()
        .collectFile ( 
            name: "prune_groups.fasta",
            storeDir: "./output/results"
        )

    //// count number of sequences passing group pruning
    ch_count_prune_groups = PRUNE_GROUPS.out.fasta.flatten().countFasta().combine(["prune_groups"])

    //// combine pruned .fasta files into a single file
    COMBINE_CHUNKS_2 ( 
        PRUNE_GROUPS.out.fasta.flatten().collect(),
        "false"
    )

    /*
    Database output
    */

    //// trim to primers based on alignment
    if ( params.trim_to_primers ) {
    
        ALIGN_OUTPUT (
            COMBINE_CHUNKS_2.out.fasta
        )

        ch_formatting_input = ALIGN_OUTPUT.out.aligned_fasta

    } else {
        ch_formatting_input = COMBINE_CHUNKS_2.out.fasta
    }

    //// format database output
    FORMAT_OUTPUT (
        ch_formatting_input,
        params.aligned_output
    )

    //// save final output
    FORMAT_OUTPUT.out.fasta
        .collectFile ( 
            name: "format_output.fasta",
            storeDir: "./output/results"
        )

    //// summarise number of taxa in database
    SUMMARISE_TAXA (
        FORMAT_OUTPUT.out.fasta
    )

    //// sequence counts for debugging
    ch_count_genbank                .view{ "${it[1]}: ${it[0]}" }
    ch_count_bold                   .view{ "${it[1]}: ${it[0]}" }
    ch_count_mito                   .view{ "${it[1]}: ${it[0]}" }
    ch_count_genome                 .view{ "${it[1]}: ${it[0]}" }
    ch_count_internal               .view{ "${it[1]}: ${it[0]}" }
    ch_count_external               .view{ "${it[1]}: ${it[0]}" }
    ch_count_input                  .view{ "${it[1]}: ${it[0]}" }
    ch_count_remove_unclassified    .view{ "${it[1]}: ${it[0]}" }
    ch_count_remove_ambiguous       .view{ "${it[1]}: ${it[0]}" }
    ch_count_filter_hmm             .view{ "${it[1]}: ${it[0]}" }
    ch_count_remove_exact           .view{ "${it[1]}: ${it[0]}" }
    ch_count_remove_tax_outliers    .view{ "${it[1]}: ${it[0]}" }
    ch_count_remove_seq_outliers    .view{ "${it[1]}: ${it[0]}" }
    ch_count_prune_groups           .view{ "${it[1]}: ${it[0]}" }

    //// collect count channels into a csv
    ch_counts_file = Channel.empty()
    ch_counts_file
        .concat ( ch_count_genbank                  .combine([1]) )
        .concat ( ch_count_bold                     .combine([2]) )
        .concat ( ch_count_mito                     .combine([3]) )
        .concat ( ch_count_genome                   .combine([4]) )
        .concat ( ch_count_internal                 .combine([5]) )
        .concat ( ch_count_external                 .combine([6]) )
        .concat ( ch_count_input                    .combine([7]) )
        .concat ( ch_count_remove_unclassified      .combine([8]) )
        .concat ( ch_count_remove_ambiguous         .combine([9]) )
        .concat ( ch_count_filter_hmm               .combine([10]) )
        .concat ( ch_count_remove_exact             .combine([11]) )
        .concat ( ch_count_remove_tax_outliers      .combine([12]) )
        .concat ( ch_count_remove_seq_outliers      .combine([13]) )
        .concat ( ch_count_prune_groups             .combine([14]) )
        .map { sequences, process, order -> "$sequences,$process,$order" }
        .collectFile ( name: 'counts.csv', seed: "sequences,process,order", newLine: true )
        .set { ch_counts_file }
    
    //// summarise the count of sequences output from stage of the pipeline
    SUMMARISE_COUNTS (
        ch_counts_file
    )

    //// train IDTAXA model
    if ( params.train_idtaxa ) {
        TRAIN_IDTAXA (
            FORMAT_OUTPUT.out.fasta
        )
    }
    
    // //// TODO: add proper model channel handling (with conditions)
    // ch_models = 
    //     Channel.empty()
    


    emit:

    // bold_db = GET_BOLD_DATABASE.out.bold_db
    ncbi_taxonomy = GET_NCBI_TAXONOMY.out.rankedlineage
    // curated_fasta = PRUNE_GROUPS.out.fasta
    // taxa_summary = TAXA_SUMMARY.out.csv
    // idtaxa_model = TRAIN_IDTAXA.out.model

}