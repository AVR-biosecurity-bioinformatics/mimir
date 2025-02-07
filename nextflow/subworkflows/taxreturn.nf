/*
 *  Run classic R::taxreturn pipeline
 */


//// modules to import
include { ADD_TAXID_GENBANK                                                 } from '../modules/add_taxid_genbank'
include { ALIGN_BATCH as ALIGN_SPECIES                                                 } from '../modules/align_batch'
include { ALIGN_CORE                                                 } from '../modules/align_core'
include { ALIGN_OTHER                                                 } from '../modules/align_other'
include { ALIGN_PRIMERS_TO_DATABASE                                                 } from '../modules/align_primers_to_database'
include { ALIGN_PRIMERS_TO_SEED                                                 } from '../modules/align_primers_to_seed'
include { ALIGN_SINGLE as ALIGN_OUTPUT                                                 } from '../modules/align_single'
include { BUILD_TRIMMED_HMM                                                 } from '../modules/build_trimmed_hmm'
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
include { GET_CORE_SEQUENCES                                         } from '../modules/get_core_sequences'
include { GET_NCBI_TAXONOMY                                         } from '../modules/get_ncbi_taxonomy'
include { HMMSEARCH_FULL                                         } from '../modules/hmmsearch_full'
include { HMMSEARCH_TRIMMED                                         } from '../modules/hmmsearch_trimmed'
include { IMPORT_INTERNAL                                                 } from '../modules/import_internal'
include { MAKE_MERGE_TABLE                                                 } from '../modules/make_merge_table'
include { MATCH_BOLD                                                 } from '../modules/match_bold'
// include { MATCH_INTERNAL                                                 } from '../modules/match_internal'
include { MERGE_ALIGNMENTS                                                 } from '../modules/merge_alignments'
include { MERGE_BOLD                                                 } from '../modules/merge_bold'
include { MERGE_SPLITS as MERGE_SPLITS_FAMILY                         } from '../modules/merge_splits'
include { MERGE_SPLITS as MERGE_SPLITS_SPECIES                             } from '../modules/merge_splits'
include { PARSE_MARKER                                                 } from '../modules/parse_marker'
include { PARSE_TARGETS                                                 } from '../modules/parse_targets'
include { PROCESS_PRIMERS                                                 } from '../modules/process_primers'
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
include { SORT_BY_LINEAGE as SORT_BY_LINEAGE_1                                                } from '../modules/sort_by_lineage'
include { SORT_BY_LINEAGE as SORT_BY_LINEAGE_2                                                } from '../modules/sort_by_lineage'
include { SPLIT_BY_RANK as SPLIT_BY_FAMILY                                                } from '../modules/split_by_rank'
include { SPLIT_BY_RANK as SPLIT_BY_SPECIES                                                } from '../modules/split_by_rank'
include { STOCKHOLM_TO_FASTA                                                 } from '../modules/stockholm_to_fasta'
include { SUMMARISE_COUNTS                                                 } from '../modules/summarise_counts'
include { SUMMARISE_TAXA                                                 } from '../modules/summarise_taxa'
include { TRAIN_IDTAXA                                                 } from '../modules/train_idtaxa'
include { TRANSLATE_SEQUENCES                                                 } from '../modules/translate_sequences'
include { TRIM_HMM_SEED                                                 } from '../modules/trim_hmm_seed'
include { TRIM_WITH_PRIMERS                                                 } from '../modules/trim_with_primers'
include { TRIM_WITH_HMM                                                 } from '../modules/trim_with_hmm'
// include { TRIM_PHMM                                                 } from '../modules/trim_phmm'



workflow TAXRETURN {

    take:

    seqs_internal


    main:

    /*
    Define channels
    */

    ch_target_taxon = channel.of ( params.target_taxon ) ?: Channel.empty()
    ch_target_rank = channel.of ( params.target_rank ) ?: Channel.of("no_rank")

    //// form ch_targets input channel
    if ( ( params.target_taxon && params.target_rank ) && !params.target_list ) { // use single taxon 
        //// combine taxon and rank into tuple for channel
        ch_target_taxon
            .combine ( ch_target_rank ) 
            .set { ch_targets }

        ch_targets.view()

    } else if ( params.target_list ) { // else use list
        //// split .csv into ch_targets
        channel.fromPath ( params.target_list, checkIfExists: true, type: 'file' )
            .splitCsv ( by: 1, header: false, strip: true )
            .map { row -> [ row[0], row[1] ] }
            .set { ch_targets }

    } else { // else invalid
        error "*** Invalid combination of '--target_taxon', '--target_rank' and '--target_list' ***"
    }

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
    // ch_phmm             = channel.fromPath( params.phmm_model, checkIfExists: true ).first()

    /*
    Set up pipeline
    */

    //// get NCBI taxonomy file
    GET_NCBI_TAXONOMY ()

    //// convert input taxon/taxa into NCBI and BOLD tax IDs
    PARSE_TARGETS (
        ch_targets,
        ch_entrez_key,
        GET_NCBI_TAXONOMY.out.synonyms,
        GET_NCBI_TAXONOMY.out.ncbi_gencodes
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

        ch_trimmed_hmm = BUILD_TRIMMED_HMM.out.hmm.first()

    } else {
        ch_trimmed_hmm = Channel.empty()
    }


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

        // //// add NCBI taxid to header of Genbank sequences 
        // ADD_TAXID_GENBANK (
        //     FETCH_GENBANK.out.fetched_seqs
        // )

        //// reformat sequence names to contain taxonomic lineage
        RENAME_GENBANK (
            FETCH_GENBANK.out.fetched_seqs,
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

        //// combine each chunk with the taxon name and rank from PARSE_TARGETS
        ch_bold_db_chunks
            .combine ( PARSE_TARGETS.out.bold ) 
            .set { ch_extract_bold_input }

        //// extract sequences from BOLD database file
        EXTRACT_BOLD (
            ch_extract_bold_input,
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
            .set { ch_translate_sequences_input }

    } else {

        ch_count_remove_unclassified = Channel.of(["NA", "remove_unclassified"])

        ch_remove_unclassified_input
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
    HMMSEARCH_FULL (
        TRANSLATE_SEQUENCES.out.translations,
        PARSE_MARKER.out.phmm
    )

    //// filter sequences by HMMSEARCH output
    FILTER_HMM (
        HMMSEARCH_FULL.out.hmmer_output,
        params.hmm_max_evalue,
        params.hmm_min_score,
        params.hmm_max_hits,
        params.hmm_min_acc,
        params.hmm_max_gap
    )

    //// combine and save intermediate file 
    if ( params.save_intermediate ) {
        FILTER_HMM.out.retained
            .map { nucleotide, protein -> nucleotide } // get just nucleotide sequences
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

    //// count number of sequences passing PHMM filter
    ch_count_filter_hmm = FILTER_HMM.out.retained.map { nucleotide, protein -> nucleotide }.countFasta().combine(["filter_hmm"]) 

    //// trim further to primer region if required
    if ( params.trim_to_primers ){
        //// do a second round of hmm searching using trimmed PHMM
        HMMSEARCH_TRIMMED (
            FILTER_HMM.out.retained,
            ch_trimmed_hmm
        )

        //// trim to trimmed HMM region
        TRIM_WITH_HMM (
            HMMSEARCH_TRIMMED.out.hmmer_output,
            params.hmm_max_evalue,
            params.hmm_min_score,
            params.hmm_max_hits,
            params.hmm_min_acc,
            params.hmm_max_gap,
            params.min_length_trimmed
        )

        //// combine and save intermediate file 
        if ( params.save_intermediate ) {
            TRIM_WITH_HMM.out.fasta
                .collectFile ( 
                    name: "trim_with_hmm.fasta",
                    storeDir: "./output/results",
                    cache: 'lenient'
                )
        }

        //// combine and save excluded sequences
        if ( params.save_intermediate ) {
            TRIM_WITH_HMM.out.removed_fasta
                .collectFile ( 
                    name: "trim_with_hmm.removed.fasta",
                    storeDir: "./output/results",
                    cache: 'lenient'
                )
        }

        //// combine and save excluded sequence table
        if ( params.save_intermediate ) {
            TRIM_WITH_HMM.out.removed_csv
                .collectFile ( 
                    name: "trim_with_hmm.removed.csv",
                    storeDir: "./output/results",
                    cache: 'lenient',
                    keepHeader: true, 
                    skip: 1
                )
        }

        //// count number of sequences passing trimming with PHMM (above min length)
        ch_count_hmm_trim = TRIM_WITH_HMM.out.fasta.countFasta().combine(["hmm_trim"]) 

        //// create output channel
        ch_hmm_output = TRIM_WITH_HMM.out.fasta 
            .filter { it.size() > 0 }
            .collect ()

    } else {
        //// null count channel
        ch_count_hmm_trim = Channel.of(["NA", "hmm_trim"])

        //// create output channel
        ch_hmm_output = FILTER_HMM.out.retained
            .map { nucleotide, protein -> nucleotide }
            .filter { it.size() > 0 } // remove empty files
            .collect ()
    }


    //// combine sequences into one .fasta file and dealign
    COMBINE_CHUNKS_1 ( 
        ch_hmm_output,
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

    //// remove ambiguous 
    if ( params.remove_ambiguous ){
        REMOVE_AMBIGUOUS (
            REMOVE_EXACT_DUPLICATES.out.fasta
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
            .set { ch_cluster_sequences_input }

    } else {

        ch_count_remove_ambiguous = Channel.of(["NA", "remove_ambiguous"])

        REMOVE_EXACT_DUPLICATES.out.fasta
            .set { ch_cluster_sequences_input }
    }

    //// cluster sequences into OTUs with mmseqs2
    CLUSTER_SEQUENCES (
        ch_cluster_sequences_input
    )

    //// remove taxonomic outliers from sequence clusters
    REMOVE_TAX_OUTLIERS (
        ch_cluster_sequences_input,
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

    //// sort .fasta by lineage string to reduce the number of files that need to be merged after splitting 
    SORT_BY_LINEAGE_1 (
        REMOVE_TAX_OUTLIERS.out.fasta
    )

    //// chunk input into SPLIT_BY_SPECIES to parallelise
    SORT_BY_LINEAGE_1.out.fasta
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
        .flatten ()
        .collect ()
        .set { ch_merge_splits_species_input }

    //// merge files from the same species across the different chunks
    MERGE_SPLITS_SPECIES (
        ch_merge_splits_species_input
    )

    //// group species-level .fasta files into batches for aligning and pruning
    ch_split_species_output.no_merge
        .flatten ()
        .mix ( MERGE_SPLITS_SPECIES.out.fasta.flatten() )
        .collect ( sort: true ) // force the channel order to be the same every time for caching -- unlikely to be a bottleneck?
        .flatten ()
        .buffer ( size: 100, remainder: true ) 
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
    if ( params.save_intermediate ) {
        PRUNE_GROUPS.out.fasta
            .flatten()
            .collectFile ( 
                name: "prune_groups.fasta",
                storeDir: "./output/results"
            )
    }

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

    // if ( params.trim_to_primers || params.aligned_output ) {
    if ( params.aligned_output ) {
        
        //// whole-database alignment method based on mafft-sparsecore.rb (adding sequences to a core alignment)

        //// define 'core' and 'other' sequences based on taxonomy and length
        GET_CORE_SEQUENCES (
            COMBINE_CHUNKS_2.out.fasta
        )

        //// align core sequences
        ALIGN_CORE (
            GET_CORE_SEQUENCES.out.core
        )

        //// add other sequences to core alignment
        ALIGN_OTHER (
            ALIGN_CORE.out.fasta,
            GET_CORE_SEQUENCES.out.other
        )

        ALIGN_OTHER.out.fasta
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

        // //// trim database to primer region
        // if ( params.trim_to_primers ){

        //     //// split into individual primer sequences for alignment
        //     PROCESS_PRIMERS.out.nuc_fasta
        //         .splitFasta ( by: 1, file: true )
        //         .set { ch_disambiguated_primers }

        //     //// add primer sequences to database alignment, output only the primer sequence (with gaps)
        //     ALIGN_PRIMERS_TO_DATABASE (
        //         ch_aligned_database,
        //         ch_disambiguated_primers
        //     )

        //     //// merge individual primer alignments into one alignment with the database
        //     ch_aligned_database
        //         .mix ( ALIGN_PRIMERS_TO_DATABASE.out.fasta )
        //         .collectFile ( name: 'merged_primers.fasta', newLine: true )
        //         .set { ch_merged_primer_alignment }

        //     //// trim alignment to primers
        //     TRIM_WITH_PRIMERS (
        //         ch_merged_primer_alignment,
        //         params.primer_fwd,
        //         params.primer_rev,
        //         params.remove_primers,
        //         params.max_primer_mismatches,
        //         params.min_length_trimmed
        //     )

        //     //// save TRIM_TO_PRIMERS output
        //     if ( params.save_intermediate ) {
        //         TRIM_WITH_PRIMERS.out.fasta
        //             .collectFile ( 
        //                 name: "trim_with_primers.fasta",
        //                 storeDir: "./output/results"
        //             )
        //     }

        //     //// save TRIM_TO_PRIMERS removed sequences
        //     if ( params.save_intermediate ) {
        //         TRIM_WITH_PRIMERS.out.removed
        //             .collectFile ( 
        //                 name: "trim_with_primers.removed.fasta",
        //                 storeDir: "./output/results"
        //             )
        //     }

        //     ch_count_trim_with_primers = TRIM_WITH_PRIMERS.out.fasta.countFasta().combine(["trim_with_primers"])

        //     ch_formatting_input = TRIM_WITH_PRIMERS.out.fasta

    } else {
        // ch_count_trim_with_primers = Channel.of(["NA", "trim_with_primers"])

        ch_formatting_input = COMBINE_CHUNKS_2.out.fasta
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
    ch_count_filter_hmm             .view{ "${it[1]}: ${it[0]}" }
    ch_count_hmm_trim               .view{ "${it[1]}: ${it[0]}" }
    ch_count_remove_exact           .view{ "${it[1]}: ${it[0]}" }
    ch_count_remove_ambiguous       .view{ "${it[1]}: ${it[0]}" }
    ch_count_remove_tax_outliers    .view{ "${it[1]}: ${it[0]}" }
    ch_count_remove_seq_outliers    .view{ "${it[1]}: ${it[0]}" }
    ch_count_prune_groups           .view{ "${it[1]}: ${it[0]}" }
    // ch_count_trim_with_primers      .view{ "${it[1]}: ${it[0]}" }

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
        .concat ( ch_count_filter_hmm               .combine([9]) )
        .concat ( ch_count_hmm_trim                 .combine([10]) )
        .concat ( ch_count_remove_exact             .combine([11]) )
        .concat ( ch_count_remove_ambiguous         .combine([12]) )
        .concat ( ch_count_remove_tax_outliers      .combine([13]) )
        .concat ( ch_count_remove_seq_outliers      .combine([14]) )
        .concat ( ch_count_prune_groups             .combine([15]) )
        // .concat ( ch_count_trim_with_primers        .combine([16]) )
        .map { sequences, process, order -> "$sequences,$process,$order" }
        .collectFile ( name: 'counts.csv', seed: "sequences,process,order", newLine: true, cache: false )
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