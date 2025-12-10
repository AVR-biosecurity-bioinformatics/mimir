/*
Filter sequences
*/


//// modules to import
include { ALIGN_BATCH as ALIGN_SPECIES                               } from '../modules/align_batch'
include { CLUSTER_SEQUENCES                                          } from '../modules/cluster_sequences'
include { COMBINE_CHUNKS                                             } from '../modules/combine_chunks'
include { FILTER_AMBIGUOUS                                           } from '../modules/filter_ambiguous'
include { FILTER_DUPLICATES                                          } from '../modules/filter_duplicates'
include { FILTER_PHMM_FULL                                           } from '../modules/filter_phmm_full'
include { FILTER_PHMM_AMPLICON                                       } from '../modules/filter_phmm_amplicon'
include { FILTER_REDUNDANT                                           } from '../modules/filter_redundant'
include { FILTER_SEQ_OUTLIERS                                        } from '../modules/filter_seq_outliers'
include { FILTER_TAX_OUTLIERS                                        } from '../modules/filter_tax_outliers'
include { FILTER_UNCLASSIFIED                                        } from '../modules/filter_unclassified'
include { HMMSEARCH_FULL                                             } from '../modules/hmmsearch_full'
include { HMMSEARCH_AMPLICON                                         } from '../modules/hmmsearch_amplicon'
include { MERGE_SPLITS as MERGE_SPLITS_GENUS                         } from '../modules/merge_splits'
include { SELECT_FINAL_SEQUENCES                                     } from '../modules/select_final_sequences'
include { SORT_BY_LINEAGE                                            } from '../modules/sort_by_lineage'
include { SPLIT_BY_RANK as SPLIT_BY_GENUS                            } from '../modules/split_by_rank'
include { TRANSLATE_SEQUENCES                                        } from '../modules/translate_sequences'


workflow FILTER_SEQUENCES {

    take:

    ch_input_seqs
    ch_internal_names
    ch_marker_coding
    ch_marker_type
    ch_gencodes
    ch_full_phmm
    ch_amplicon_phmm
    ch_primer_info

    main:
      
    /*
    Sequence filtering
    */

    //// remove empty .fasta files
    ch_input_seqs
        .filter { it.size() > 0 }
        .set { ch_translate_sequences_input }

    //// translate sequences in all six frames
    TRANSLATE_SEQUENCES (
        ch_translate_sequences_input,
        ch_marker_coding,
        ch_marker_type,
        ch_gencodes
    )

    //// search translated sequences for hits against PHMM model
    HMMSEARCH_FULL (
        TRANSLATE_SEQUENCES.out.translations,
        ch_full_phmm
    )

    //// filter sequences by HMMSEARCH output
    FILTER_PHMM_FULL (
        HMMSEARCH_FULL.out.hmmer_output,
        params.phmm_max_evalue,
        params.phmm_min_score,
        params.phmm_max_hits,
        params.phmm_min_acc,
        params.phmm_max_gap,
        params.phmm_min_length,
        params.phmm_min_cov
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


    //// trim and/or filter based on amplicon PHMM
    if ( ch_amplicon_phmm ){
        //// do a second round of hmm searching using trimmed PHMM
        HMMSEARCH_AMPLICON (
            FILTER_PHMM_FULL.out.retained,
            ch_amplicon_phmm
        )

        //// apply HMMER hits to nucleotide sequences
        FILTER_PHMM_AMPLICON (
            HMMSEARCH_AMPLICON.out.hmmer_output,
            ch_primer_info.first(),
            params.trim_to_amplicon,
            params.amplicon_min_length,
            params.amplicon_min_cov,
            params.remove_primers
        )

        //// combine and save intermediate file 
        if ( params.save_intermediate ) {
            FILTER_PHMM_AMPLICON.out.fasta
                .collectFile ( 
                    name: "filter_phmm_amplicon.fasta",
                    storeDir: "./output/results",
                    cache: 'lenient'
                )
        }

        //// combine and save excluded sequences
        if ( params.save_intermediate ) {
            FILTER_PHMM_AMPLICON.out.removed_fasta
                .collectFile ( 
                    name: "filter_phmm_amplicon.removed.fasta",
                    storeDir: "./output/results",
                    cache: 'lenient'
                )
        }

        //// combine and save excluded sequence table
        if ( params.save_intermediate ) {
            FILTER_PHMM_AMPLICON.out.removed_csv
                .collectFile ( 
                    name: "filter_phmm_amplicon.removed.csv",
                    storeDir: "./output/results",
                    cache: 'lenient',
                    keepHeader: true, 
                    skip: 1
                )
        }

        //// count number of sequences passing trimming with PHMM (above min length)
        ch_count_filter_phmm_amplicon = FILTER_PHMM_AMPLICON.out.fasta.countFasta().combine(["filter_phmm_amplicon"]) 

        //// sequence names that failed filter
        FILTER_PHMM_AMPLICON.out.removed_fasta
            .collectFile ( name: 'filter_phmm_amplicon.fasta', newLine: true, cache: false )
            .set { ch_fates_filter_phmm_amplicon }

        //// create output channel
        ch_hmm_output = FILTER_PHMM_AMPLICON.out.fasta 
            .filter { it.size() > 0 }
            .collect ()

    } else {
        //// null count channel
        ch_count_filter_phmm_amplicon = Channel.of(["NA", "filter_phmm_amplicon"])

        //// empty fate channek
        ch_fates_filter_phmm_amplicon = Channel.empty()

        //// create output channel
        ch_hmm_output = FILTER_PHMM_FULL.out.retained
            .map { nucleotide, protein -> nucleotide }
            .filter { it.size() > 0 } // remove empty files
            .collect ()
    }


    //// combine sequences into one .fasta file and dealign
    COMBINE_CHUNKS ( 
        ch_hmm_output,
        "true"
    )
 
    //// remove exact duplicate sequences if they exist
    FILTER_DUPLICATES (
        COMBINE_CHUNKS.out.fasta
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
            .set { ch_sort_input }

    } else {

        ch_count_filter_ambiguous = Channel.of(["NA", "filter_ambiguous"])

        ch_fates_filter_ambiguous = Channel.empty()

        FILTER_DUPLICATES.out.fasta
            .set { ch_sort_input }
    }

    
    
    //////////////////////////////////////////////// main changes start here

    //// sort .fasta by lineage string to reduce the number of files that need to be merged after splitting 
    SORT_BY_LINEAGE (
        ch_sort_input
    )

    //// chunk input into SPLIT_BY_GENUS to parallelise
    SORT_BY_LINEAGE.out.fasta
        .splitFasta ( by: params.split_rank_chunk, file: true )
        .set { ch_split_genus_input }

    //// split .fasta by taxonomic lineage down to genus level
    SPLIT_BY_GENUS (
        ch_split_genus_input,
        "genus"
    )

    //// group files with the same name
    SPLIT_BY_GENUS.out.fasta
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
        .set { ch_split_genus_output }

    //// group merging files into groups of 1000 genera for merging
    ch_split_genus_output.merge
        .buffer ( size: 1000, remainder: true )
        .flatten ()
        .collect ()
        .set { ch_merge_splits_genus_input }

    //// merge files from the same species across the different chunks
    MERGE_SPLITS_GENUS (
        ch_merge_splits_genus_input
    )

    //// group genus-level .fasta files into batches for redundancy filtering
    ch_split_genus_output.no_merge
        .flatten ()
        .mix ( MERGE_SPLITS_GENUS.out.fasta.flatten() )
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





    // //// cluster sequences into OTUs with mmseqs2
    // CLUSTER_SEQUENCES (
    //     ch_cluster_sequences_input,
    //     params.cluster_threshold
    // )

    // //// remove taxonomic outliers from sequence clusters
    // FILTER_TAX_OUTLIERS (
    //     ch_cluster_sequences_input,
    //     CLUSTER_SEQUENCES.out.tsv,
    //     params.cluster_rank,
    //     params.cluster_threshold,
    //     params.cluster_confidence
    // )

    // //// combine and save intermediate file 
    // if ( params.save_intermediate ) {
    //     FILTER_TAX_OUTLIERS.out.fasta
    //         .collectFile ( 
    //             name: "filter_tax_outliers.fasta",
    //             storeDir: "./output/results"
    //         )
    // }

    // //// count number of sequences passing taxonomic decontamination
    // ch_count_filter_tax_outliers = FILTER_TAX_OUTLIERS.out.fasta.countFasta().combine(["filter_tax_outliers"])

    // //// sequence names that failed filter
    // FILTER_TAX_OUTLIERS.out.removed
    //     .collectFile ( name: 'filter_tax_outliers.fasta', newLine: true, cache: false )
    //     .set { ch_fates_filter_tax_outliers }

    

    // //// align species-level .fasta in batches
    // ALIGN_SPECIES (
    //     FILTER_REDUNDANT.out.fasta
    // )

    // //// remove sequence outliers from species clusters
    // FILTER_SEQ_OUTLIERS (
    //     ALIGN_SPECIES.out.aligned_fasta,
    //     params.dist_threshold
    // )

    // //// combine and save intermediate file 
    // if ( params.save_intermediate ) {
    //     FILTER_SEQ_OUTLIERS.out.retained_fasta
    //         .flatten()
    //         .collectFile ( 
    //             name: "filter_seq_outliers.fasta",
    //             storeDir: "./output/results"
    //         )
    // }

    // //// count number of sequences passing taxonomic decontamination
    // ch_count_filter_seq_outliers = FILTER_SEQ_OUTLIERS.out.retained_fasta.flatten().countFasta().combine(["filter_seq_outliers"])

    // //// sequence names that failed filter
    // FILTER_SEQ_OUTLIERS.out.removed
    //     .flatten()
    //     .collectFile ( name: 'filter_seq_outliers.fasta', newLine: true, cache: false )
    //     .set { ch_fates_filter_seq_outliers }

    // //// remove unclassified sequences
    // if ( params.remove_unclassified == "all_ranks" || params.remove_unclassified == "any_ranks" || params.remove_unclassified == "terminal" ) {
    //     FILTER_UNCLASSIFIED (
    //         FILTER_SEQ_OUTLIERS.out.retained_fasta,
    //         params.remove_unclassified
    //     )        

    //     ch_count_filter_unclassified = FILTER_UNCLASSIFIED.out.fasta.countFasta().combine(["filter_unclassified"])

    //     //// combine and save intermediate file 
    //     if ( params.save_intermediate ) {
    //         FILTER_UNCLASSIFIED.out.fasta
    //             .collectFile ( 
    //                 name: "filter_unclassified.fasta",
    //                 storeDir: "./output/results",
    //                 cache: 'lenient'
    //             )
    //     }

    //     //// sequence names that failed filter
    //     FILTER_UNCLASSIFIED.out.removed
    //         .collectFile ( name: 'filter_unclassified.fasta', newLine: true, cache: false )
    //         .set { ch_fates_filter_unclassified }

    //     FILTER_UNCLASSIFIED.out.fasta
    //         .filter { it.size() > 0 }
    //         .set { ch_select_final_input }

    // } else {

    //     ch_count_filter_unclassified = Channel.of(["NA", "filter_unclassified"])

    //     ch_fates_filter_unclassified = Channel.empty()

    //     FILTER_SEQ_OUTLIERS.out.retained_fasta
    //         .filter { it.size() > 0 }
    //         .set { ch_select_final_input }
    
    // }

    // //// remove reundant sequences from each species, preferentially retaining internal sequences
    // SELECT_FINAL_SEQUENCES (
    //     ch_select_final_input,
    //     ch_internal_names,
    //     params.max_group_size,
    //     params.selection_method
    // )

    // //// save SELECT_FINAL_SEQUENCES output
    // if ( params.save_intermediate ) {
    //     SELECT_FINAL_SEQUENCES.out.fasta
    //         .flatten()
    //         .collectFile ( 
    //             name: "select_final_sequences.fasta",
    //             storeDir: "./output/results"
    //         )
    // }

    // //// count number of sequences passing group pruning
    // ch_count_select_final_sequences = SELECT_FINAL_SEQUENCES.out.fasta.flatten().countFasta().combine(["select_final_sequences"])

    // //// sequence names that failed filter
    // SELECT_FINAL_SEQUENCES.out.removed
    //     .collectFile ( name: 'select_final_sequences.fasta', newLine: true, cache: false )
    //     .set { ch_fates_select_final_sequences }

    // //// sequence names included in final datasase
    // SELECT_FINAL_SEQUENCES.out.fasta
    //     .collectFile ( name: 'final_database.fasta', newLine: true, cache: false )
    //     .set { ch_fates_final_database }

    
    // /// collect sequence fates into a list of .csv files for concatenation
    // Channel.empty()
    //     .concat ( ch_fates_filter_unclassified )
    //     .concat ( ch_fates_filter_phmm_full )
    //     .concat ( ch_fates_filter_phmm_amplicon )
    //     .concat ( ch_fates_filter_duplicates )
    //     .concat ( ch_fates_filter_ambiguous )
    //     .concat ( ch_fates_filter_tax_outliers )
    //     .concat ( ch_fates_filter_redundant )
    //     .concat ( ch_fates_filter_seq_outliers )
    //     .concat ( ch_fates_select_final_sequences )
    //     .concat ( ch_fates_final_database )
    //     .collect ()
    //     .set { ch_fates }

    emit:

    // ch_final_seqs = SELECT_FINAL_SEQUENCES.out.fasta
    ch_final_seqs = Channel.empty()
    ch_fates = Channel.empty()
    ch_count_filter_unclassified = Channel.empty()
    ch_count_filter_phmm_full = Channel.empty()
    ch_count_filter_phmm_amplicon = Channel.empty()
    ch_count_filter_duplicates = Channel.empty()
    ch_count_filter_ambiguous = Channel.empty()
    ch_count_filter_tax_outliers = Channel.empty()
    ch_count_filter_redundant = Channel.empty()
    ch_count_filter_seq_outliers = Channel.empty()
    ch_count_select_final_sequences = Channel.empty()

}