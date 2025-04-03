/*
 *  Run classic R::taxreturn pipeline
 */


//// modules to import
include { ALIGN_BATCH as ALIGN_SPECIES                               } from '../modules/align_batch'
include { ALIGN_CORE                                                 } from '../modules/align_core'
include { ALIGN_OTHER                                                } from '../modules/align_other'
include { ALIGN_PRIMERS_TO_SEED                                      } from '../modules/align_primers_to_seed'
include { BUILD_TRIMMED_HMM                                          } from '../modules/build_trimmed_hmm'
include { CHECK_KEY_SPECIES                                          } from '../modules/check_key_species'
include { CLUSTER_SEQUENCES                                          } from '../modules/cluster_sequences'
include { COMBINE_CHUNKS as COMBINE_CHUNKS_1                         } from '../modules/combine_chunks'
include { COMBINE_CHUNKS as COMBINE_CHUNKS_2                         } from '../modules/combine_chunks'
include { CONCAT_CSV as CONCAT_FATES                                 } from '../modules/concat_csv'
include { CONCAT_CSV as CONCAT_SOURCES                               } from '../modules/concat_csv'
include { EXTRACT_BOLD                                               } from '../modules/extract_bold'
include { FETCH_GENBANK                                              } from '../modules/fetch_genbank'
include { FILTER_AMBIGUOUS                                           } from '../modules/filter_ambiguous'
include { FILTER_DUPLICATES                                          } from '../modules/filter_duplicates'
include { FILTER_PHMM_FULL                                           } from '../modules/filter_phmm_full'
include { FILTER_PHMM_TRIMMED                                        } from '../modules/filter_phmm_trimmed'
include { FILTER_REDUNDANT                                           } from '../modules/filter_redundant'
include { FILTER_SEQ_OUTLIERS                                        } from '../modules/filter_seq_outliers'
include { FILTER_TAX_OUTLIERS                                        } from '../modules/filter_tax_outliers'
include { FILTER_UNCLASSIFIED                                        } from '../modules/filter_unclassified'
include { FORMAT_OUTPUT                                              } from '../modules/format_output'
include { GET_BOLD_DATABASE                                          } from '../modules/get_bold_database'
include { GET_CORE_SEQUENCES                                         } from '../modules/get_core_sequences'
include { GET_NCBI_TAXONOMY                                          } from '../modules/get_ncbi_taxonomy'
include { HMMSEARCH_FULL                                             } from '../modules/hmmsearch_full'
include { HMMSEARCH_TRIMMED                                          } from '../modules/hmmsearch_trimmed'
include { IMPORT_INTERNAL                                            } from '../modules/import_internal'
include { JOIN_SOURCES_FATES                                         } from '../modules/join_sources_fates'
include { MATCH_BOLD                                                 } from '../modules/match_bold'
include { MERGE_BOLD                                                 } from '../modules/merge_bold'
include { MERGE_SPLITS as MERGE_SPLITS_SPECIES                       } from '../modules/merge_splits'
include { PARSE_MARKER                                               } from '../modules/parse_marker'
include { PARSE_TARGETS                                              } from '../modules/parse_targets'
include { PROCESS_PRIMERS                                            } from '../modules/process_primers'
include { QUERY_GENBANK                                              } from '../modules/query_genbank'
include { RENAME_GENBANK                                             } from '../modules/rename_genbank'
include { SELECT_FINAL_SEQUENCES                                     } from '../modules/select_final_sequences'
include { SEQUENCE_TRACKER                                           } from '../modules/sequence_tracker'
include { SORT_BY_LINEAGE as SORT_BY_LINEAGE_1                       } from '../modules/sort_by_lineage'
include { SPLIT_BY_RANK as SPLIT_BY_SPECIES                          } from '../modules/split_by_rank'
include { STOCKHOLM_TO_FASTA                                         } from '../modules/stockholm_to_fasta'
include { TRAIN_IDTAXA                                               } from '../modules/train_idtaxa'
include { TRANSLATE_SEQUENCES                                        } from '../modules/translate_sequences'
include { TRIM_HMM_SEED                                              } from '../modules/trim_hmm_seed'
include { VALIDATE_KEY_SPECIES                                       } from '../modules/validate_key_species'



workflow TAXRETURN {

    take:

    seqs_internal


    main:

    /*
    Define channels
    */

    ch_target_taxon = channel.of ( params.target_taxon ) ?: Channel.empty()
    ch_target_rank = channel.of ( params.target_rank ) ?: Channel.of("no_rank")
    
    if ( params.key_species_list ){
        ch_key_species_list = channel.fromPath ( params.key_species_list, checkIfExists: true, type: 'file' )
    } else {
        ch_key_species_list = Channel.empty()
    }

    //// form ch_targets input channel
    if ( ( params.target_taxon && params.target_rank ) && !params.target_list ) { // use single taxon 
        //// combine taxon and rank into tuple for channel
        ch_target_taxon
            .combine ( ch_target_rank ) 
            .set { ch_targets }

    } else if ( params.target_list ) { // else use list
        //// split .csv into ch_targets
        channel.fromPath ( params.target_list, checkIfExists: true, type: 'file' )
            .splitCsv ( by: 1, header: false, strip: true )
            .unique()
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
    GET_NCBI_TAXONOMY (
        "dummy"
    )

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

        //// reformat sequence names to contain taxonomic lineage
        RENAME_GENBANK (
            FETCH_GENBANK.out.fetched_seqs,
            GET_NCBI_TAXONOMY.out.rankedlineage_noname,
            params.placeholder_as_unclassified
        )

        //// populate empty ch_genbank_fasta channel
        RENAME_GENBANK.out.fasta
            .set { ch_genbank_fasta }
    }
    

    //// count number of sequences downloaded from Genbank
    ch_count_genbank = ch_genbank_fasta.countFasta().combine(["genbank"])

    //// collect genbank sequences
    ch_genbank_fasta
        .collectFile ( name: 'genbank.fasta', newLine: true, cache: true )
        .set { ch_source_genbank }


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
            GET_NCBI_TAXONOMY.out.synonyms,
            params.placeholder_as_unclassified
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

    } 

    //// count number of BOLD sequences
    ch_count_bold = ch_bold_fasta.countFasta().combine(["bold"])

    //// save names of BOLD sequences
    ch_bold_fasta
        .collectFile ( name: 'bold.fasta', newLine: true, cache: true )
        .set { ch_source_bold }

    //// count number of mitochondrial sequences 
    ch_count_mito = ch_mito_fasta.countFasta().combine(["mito"])

    //// save names of mito sequences
    ch_mito_fasta
        .collectFile ( name: 'mito.fasta', newLine: true, cache: true )
        .set { ch_source_mito }

    //// count number of genome assembly-derived sequences 
    ch_count_genome = ch_genome_fasta.countFasta().combine(["genome"])

    //// save names of genome sequences
    ch_genome_fasta
        .collectFile ( name: 'genome.fasta', newLine: true, cache: true )
        .set { ch_source_genome }



    /*
    Import and format internal sequences
    */

    if ( ch_internal_seqs ){
        //// import internal sequences from file and check format
        IMPORT_INTERNAL (
            ch_internal_seqs
        )
        
        //// populate and chunk internal channel
        IMPORT_INTERNAL.out.fasta
            .splitText( by: params.input_chunk_size, file: true )
            .set { ch_internal_fasta }

        //// save names of internal sequences
        ch_internal_fasta
            .collectFile ( name: 'internal.fasta', newLine: true, cache: true )
            .set { ch_source_internal }

    } else {
        ch_source_internal = Channel.empty()
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

    // //// save names of all input sequences
    // ch_input_seqs
    //     .flatten()
    //     .splitFasta ( by: 1, record: [ header:true ] )
    //     .map { record -> record.header }
    //     .collectFile ( name: 'input_names.txt', newLine: true, cache: true )
    //     .set { ch_input_names_file }


    
    /*
    Sequence filtering
    */

    //// remove unclassified sequences
    if ( params.remove_unclassified == "all_ranks" || params.remove_unclassified == "any_ranks" || params.remove_unclassified == "terminal" ) {
        FILTER_UNCLASSIFIED (
            ch_input_seqs,
            params.remove_unclassified
        )        

        ch_count_filter_unclassified = FILTER_UNCLASSIFIED.out.fasta.countFasta().combine(["filter_unclassified"])

        //// combine and save intermediate file 
        if ( params.save_intermediate ) {
            FILTER_UNCLASSIFIED.out.fasta
                .collectFile ( 
                    name: "filter_unclassified.fasta",
                    storeDir: "./output/results",
                    cache: 'lenient'
                )
        }

        //// sequence names that failed filter
        FILTER_UNCLASSIFIED.out.removed
            .collectFile ( name: 'filter_unclassified.fasta', newLine: true, cache: false )
            .set { ch_fates_filter_unclassified }

        FILTER_UNCLASSIFIED.out.fasta
            .filter { it.size() > 0 }
            .set { ch_translate_sequences_input }

    } else {

        ch_count_filter_unclassified = Channel.of(["NA", "filter_unclassified"])

        ch_fates_filter_unclassified = Channel.empty()

        ch_input_seqs
            .filter { it.size() > 0 }
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
    FILTER_PHMM_FULL (
        HMMSEARCH_FULL.out.hmmer_output,
        params.hmm_max_evalue,
        params.hmm_min_score,
        params.hmm_max_hits,
        params.hmm_min_acc,
        params.hmm_max_gap
    )

    //// combine and save intermediate file 
    if ( params.save_intermediate ) {
        FILTER_PHMM_FULL.out.retained
            .map { nucleotide, protein -> nucleotide } // get just nucleotide sequences
            .collectFile ( 
                name: "filter_phmm_full.fasta",
                storeDir: "./output/results",
                cache: 'lenient'
            )
    }

    //// combine and save excluded sequences
    if ( params.save_intermediate ) {
        FILTER_PHMM_FULL.out.removed_fasta
            .collectFile ( 
                name: "filter_phmm_full.removed.fasta",
                storeDir: "./output/results",
                cache: 'lenient'
            )
    }

    //// combine and save excluded sequence table
    if ( params.save_intermediate ) {
        FILTER_PHMM_FULL.out.removed_csv
            .collectFile ( 
                name: "filter_phmm_full.removed.csv",
                storeDir: "./output/results",
                cache: 'lenient',
                keepHeader: true, 
                skip: 1
            )
    }

    //// count number of sequences passing PHMM filter
    ch_count_filter_phmm_full = FILTER_PHMM_FULL.out.retained.map { nucleotide, protein -> nucleotide }.countFasta().combine(["filter_phmm_full"]) 

    //// sequence names that failed filter
    FILTER_PHMM_FULL.out.removed_fasta
        .collectFile ( name: 'filter_phmm_full.fasta', newLine: true, cache: false )
        .set { ch_fates_filter_phmm_full }


    //// trim further to primer region if required
    if ( params.trim_to_primers ){
        //// do a second round of hmm searching using trimmed PHMM
        HMMSEARCH_TRIMMED (
            FILTER_PHMM_FULL.out.retained,
            ch_trimmed_hmm
        )

        //// trim to trimmed HMM region
        FILTER_PHMM_TRIMMED (
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
            FILTER_PHMM_TRIMMED.out.fasta
                .collectFile ( 
                    name: "filter_phmm_trimmed.fasta",
                    storeDir: "./output/results",
                    cache: 'lenient'
                )
        }

        //// combine and save excluded sequences
        if ( params.save_intermediate ) {
            FILTER_PHMM_TRIMMED.out.removed_fasta
                .collectFile ( 
                    name: "filter_phmm_trimmed.removed.fasta",
                    storeDir: "./output/results",
                    cache: 'lenient'
                )
        }

        //// combine and save excluded sequence table
        if ( params.save_intermediate ) {
            FILTER_PHMM_TRIMMED.out.removed_csv
                .collectFile ( 
                    name: "filter_phmm_trimmed.removed.csv",
                    storeDir: "./output/results",
                    cache: 'lenient',
                    keepHeader: true, 
                    skip: 1
                )
        }

        //// count number of sequences passing trimming with PHMM (above min length)
        ch_count_filter_phmm_trimmed = FILTER_PHMM_TRIMMED.out.fasta.countFasta().combine(["filter_phmm_trimmed"]) 

        //// sequence names that failed filter
        FILTER_PHMM_TRIMMED.out.removed_fasta
            .collectFile ( name: 'filter_phmm_trimmed.fasta', newLine: true, cache: false )
            .set { ch_fates_filter_phmm_trimmed }

        //// create output channel
        ch_hmm_output = FILTER_PHMM_TRIMMED.out.fasta 
            .filter { it.size() > 0 }
            .collect ()

    } else {
        //// null count channel
        ch_count_filter_phmm_trimmed = Channel.of(["NA", "filter_phmm_trimmed"])

        //// empty fate channek
        ch_fates_filter_phmm_trimmed = Channel.empty()

        //// create output channel
        ch_hmm_output = FILTER_PHMM_FULL.out.retained
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
    FILTER_DUPLICATES (
        COMBINE_CHUNKS_1.out.fasta
    )

    //// combine and save intermediate file 
    if ( params.save_intermediate ) {
        FILTER_DUPLICATES.out.fasta
            .collectFile ( 
                name: "remove_exact_duplicates.fasta",
                storeDir: "./output/results"
            )
    }

    //// count number of sequences passing exact deduplication
    ch_count_filter_duplicates = FILTER_DUPLICATES.out.fasta.countFasta().combine(["filter_duplicates"])

    //// sequence names that failed filter
    FILTER_DUPLICATES.out.removed
        .collectFile ( name: 'filter_duplicates.fasta', newLine: true, cache: false )
        .set { ch_fates_filter_duplicates }

    //// remove ambiguous 
    if ( params.remove_ambiguous ){
        FILTER_AMBIGUOUS (
            FILTER_DUPLICATES.out.fasta
        )

        ch_count_filter_ambiguous = FILTER_AMBIGUOUS.out.fasta.countFasta().combine(["filter_ambiguous"])

        //// combine and save intermediate file 
        if ( params.save_intermediate ) {
            FILTER_AMBIGUOUS.out.fasta
                .collectFile ( 
                    name: "filter_ambiguous.fasta",
                    storeDir: "./output/results",
                    cache: 'lenient'
                )
        }

        //// sequence names that failed filter
        FILTER_AMBIGUOUS.out.removed
            .collectFile ( name: 'filter_ambiguous.fasta', newLine: true, cache: false )
            .set { ch_fates_filter_ambiguous }

        FILTER_AMBIGUOUS.out.fasta
            .set { ch_cluster_sequences_input }

    } else {

        ch_count_filter_ambiguous = Channel.of(["NA", "filter_ambiguous"])

        ch_fates_filter_ambiguous = Channel.empty()

        FILTER_DUPLICATES.out.fasta
            .set { ch_cluster_sequences_input }
    }

    //// cluster sequences into OTUs with mmseqs2
    CLUSTER_SEQUENCES (
        ch_cluster_sequences_input
    )

    //// remove taxonomic outliers from sequence clusters
    FILTER_TAX_OUTLIERS (
        ch_cluster_sequences_input,
        CLUSTER_SEQUENCES.out.tsv
    )

    //// combine and save intermediate file 
    if ( params.save_intermediate ) {
        FILTER_TAX_OUTLIERS.out.fasta
            .collectFile ( 
                name: "filter_tax_outliers.fasta",
                storeDir: "./output/results"
            )
    }

    //// count number of sequences passing taxonomic decontamination
    ch_count_filter_tax_outliers = FILTER_TAX_OUTLIERS.out.fasta.countFasta().combine(["filter_tax_outliers"])

    //// sequence names that failed filter
    FILTER_TAX_OUTLIERS.out.removed
        .collectFile ( name: 'filter_tax_outliers.fasta', newLine: true, cache: false )
        .set { ch_fates_filter_tax_outliers }

    //// sort .fasta by lineage string to reduce the number of files that need to be merged after splitting 
    SORT_BY_LINEAGE_1 (
        FILTER_TAX_OUTLIERS.out.fasta
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
        .buffer ( size: 500, remainder: true ) 
        .set { ch_filter_redundant_input }

    //// filter out redundant (ie. identical and contained) sequences within each species, counting the number of sequences absorbed
    FILTER_REDUNDANT (
        ch_filter_redundant_input
    )

    //// save FILTER_REDUNDANT output
    if ( params.save_intermediate ) {
        FILTER_REDUNDANT.out.fasta.map{fasta, csv -> fasta}
            .flatten()
            .collectFile ( 
                name: "filter_redundant.fasta",
                storeDir: "./output/results"
            )
    }

    //// count number of sequences passing group pruning
    ch_count_filter_redundant = FILTER_REDUNDANT.out.fasta.map{fasta, csv -> fasta}.flatten().countFasta().combine(["filter_redundant"])

    //// sequence names that failed filter
    FILTER_REDUNDANT.out.removed
        .flatten()
        .collectFile ( name: 'filter_redundant.fasta', newLine: true, cache: false )
        .set { ch_fates_filter_redundant }

    //// align species-level .fasta in batches
    ALIGN_SPECIES (
        FILTER_REDUNDANT.out.fasta
    )

    //// remove sequence outliers from species clusters
    FILTER_SEQ_OUTLIERS (
        ALIGN_SPECIES.out.aligned_fasta,
        params.dist_threshold
    )

    //// combine and save intermediate file 
    if ( params.save_intermediate ) {
        FILTER_SEQ_OUTLIERS.out.retained_fasta
            .flatten()
            .collectFile ( 
                name: "filter_seq_outliers.fasta",
                storeDir: "./output/results"
            )
    }

    //// count number of sequences passing taxonomic decontamination
    ch_count_filter_seq_outliers = FILTER_SEQ_OUTLIERS.out.retained_fasta.flatten().countFasta().combine(["filter_seq_outliers"])

    //// sequence names that failed filter
    FILTER_SEQ_OUTLIERS.out.removed
        .flatten()
        .collectFile ( name: 'filter_seq_outliers.fasta', newLine: true, cache: false )
        .set { ch_fates_filter_seq_outliers }

    //// remove reundant sequences from each species, preferentially retaining internal sequences
    SELECT_FINAL_SEQUENCES (
        FILTER_SEQ_OUTLIERS.out.retained_fasta,
        ch_internal_names,
        params.max_group_size
    )

    //// save SELECT_FINAL_SEQUENCES output
    if ( params.save_intermediate ) {
        SELECT_FINAL_SEQUENCES.out.fasta
            .flatten()
            .collectFile ( 
                name: "select_final_sequences.fasta",
                storeDir: "./output/results"
            )
    }

    //// count number of sequences passing group pruning
    ch_count_select_final_sequences = SELECT_FINAL_SEQUENCES.out.fasta.flatten().countFasta().combine(["select_final_sequences"])

    //// sequence names that failed filter
    SELECT_FINAL_SEQUENCES.out.removed
        .collectFile ( name: 'select_final_sequences.fasta', newLine: true, cache: false )
        .set { ch_fates_select_final_sequences }

    //// sequence names included in final datasase
    SELECT_FINAL_SEQUENCES.out.fasta
        .collectFile ( name: 'final_database.fasta', newLine: true, cache: false )
        .set { ch_fates_final_database }

    //// combine pruned .fasta files into a single file
    COMBINE_CHUNKS_2 ( 
        SELECT_FINAL_SEQUENCES.out.fasta.flatten().collect(),
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

    } else {

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

    //// sequence counts for debugging
    if ( params.debug_mode ){
        ch_count_genbank                                .view{ "${it[1]}: ${it[0]}" }
        ch_count_bold                                   .view{ "${it[1]}: ${it[0]}" }
        ch_count_mito                                   .view{ "${it[1]}: ${it[0]}" }
        ch_count_genome                                 .view{ "${it[1]}: ${it[0]}" }
        ch_count_internal                               .view{ "${it[1]}: ${it[0]}" }
        ch_count_external                               .view{ "${it[1]}: ${it[0]}" }
        ch_count_input                                  .view{ "${it[1]}: ${it[0]}" }
        ch_count_filter_unclassified                    .view{ "${it[1]}: ${it[0]}" }
        ch_count_filter_phmm_full                       .view{ "${it[1]}: ${it[0]}" }
        ch_count_filter_phmm_trimmed                    .view{ "${it[1]}: ${it[0]}" }
        ch_count_filter_duplicates                      .view{ "${it[1]}: ${it[0]}" }
        ch_count_filter_ambiguous                       .view{ "${it[1]}: ${it[0]}" }
        ch_count_filter_tax_outliers                    .view{ "${it[1]}: ${it[0]}" }
        ch_count_filter_redundant                       .view{ "${it[1]}: ${it[0]}" }
        ch_count_filter_seq_outliers                    .view{ "${it[1]}: ${it[0]}" }
        ch_count_select_final_sequences                       .view{ "${it[1]}: ${it[0]}" }
    }

    //// collect sequence origins into a list of .csv files for concatenation
    Channel.empty()
        .concat ( ch_source_genbank )
        .concat ( ch_source_bold )
        .concat ( ch_source_mito )
        .concat ( ch_source_genome )
        .concat ( ch_source_internal )
        .collect ()
        .set { ch_sources }

    //// concat sequence sources into one .csv
    CONCAT_SOURCES (
        ch_sources,
        "sources"
    )

    //// collect sequence fates into a list of .csv files for concatenation
    Channel.empty()
        .concat ( ch_fates_filter_unclassified )
        .concat ( ch_fates_filter_phmm_full )
        .concat ( ch_fates_filter_phmm_trimmed )
        .concat ( ch_fates_filter_duplicates )
        .concat ( ch_fates_filter_ambiguous )
        .concat ( ch_fates_filter_tax_outliers )
        .concat ( ch_fates_filter_redundant )
        .concat ( ch_fates_filter_seq_outliers )
        .concat ( ch_fates_select_final_sequences )
        .concat ( ch_fates_final_database )
        .collect ()
        .set { ch_fates }

    //// concat sequence fates into one .csv
    CONCAT_FATES (
        ch_fates,
        "fates"
    )

    //// join fates and sources together into a single .csv
    JOIN_SOURCES_FATES (
        CONCAT_SOURCES.out.csv,
        CONCAT_FATES.out.csv
    )

    //// process the source and fate of every sequences through the pipeline
    SEQUENCE_TRACKER (
        JOIN_SOURCES_FATES.out.csv
    )


    if ( params.key_species_list ){

        //// check for sequences belonging to key species
        CHECK_KEY_SPECIES (
            SEQUENCE_TRACKER.out.sf_meta,
            ch_key_species_list
        )

        ch_key_sequence_names = CHECK_KEY_SPECIES.out.seq_names.flatten()

        //// validate sequences from 
        VALIDATE_KEY_SPECIES (
            ch_key_sequence_names,
            FORMAT_OUTPUT.out.fasta,
            params.add_root
        )
        
    }

    //// train IDTAXA model
    if ( params.train_idtaxa ) {
        TRAIN_IDTAXA (
            FORMAT_OUTPUT.out.fasta,
            params.idtaxa_max_group_size,
            params.idtaxa_max_iterations,
            params.idtaxa_allow_group_removal
        )
    }
    


    emit:

    // bold_db = GET_BOLD_DATABASE.out.bold_db
    ncbi_taxonomy = GET_NCBI_TAXONOMY.out.rankedlineage
    // curated_fasta = PRUNE_GROUPS.out.fasta
    // taxa_summary = TAXA_SUMMARY.out.csv
    // idtaxa_model = TRAIN_IDTAXA.out.model

}