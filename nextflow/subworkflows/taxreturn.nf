/*
 *  Run classic R::taxreturn pipeline
 */


//// modules to import
include { ADD_TAXID_GENBANK                                                 } from '../modules/add_taxid_genbank'
include { COMBINE_CHUNKS                                                 } from '../modules/combine_chunks'
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
include { QUERY_GENBANK                                                 } from '../modules/query_genbank'
include { REFORMAT_NAMES                                                 } from '../modules/reformat_names'
include { REMOVE_CONTAM                                                 } from '../modules/remove_contam'
include { REMOVE_EXACT_DUPLICATES                                                 } from '../modules/remove_exact_duplicates'
include { RENAME_GENBANK                                                 } from '../modules/rename_genbank'
include { RESOLVE_SYNONYMS                                                 } from '../modules/resolve_synonyms'
include { SUMMARISE_COUNTS                                                 } from '../modules/summarise_counts'
include { SUMMARISE_TAXA                                                 } from '../modules/summarise_taxa'
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
    ch_genbank_fasta          = Channel.empty()
    ch_bold_fasta             = Channel.empty()
    ch_mito_fasta             = Channel.empty()
    ch_genome_fasta           = Channel.empty()
    ch_internal_fasta         = Channel.empty()

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
    PARSE_TARGETS.out.taxon_id_rank
        .map { id, rank -> 
            [ id.text.trim(), rank ] }
        .set { ch_taxon_idrank }

    // //// parse marker gene into formats understandable for each database
    // PARSE_MARKER (

    // )

    /*
    Getting Genbank sequences
    */

    //// query GenBank to get list of nucleotide IDs (including mitochondrial genomes if requested)
    QUERY_GENBANK (
        ch_taxon_idrank,
        "COI[GENE] OR COX1[GENE] OR COXI[GENE]"
    )

    //// chunk list of accessions for fetching
    QUERY_GENBANK.out.seq_acc
        .splitText( by: params.chunk_size, file: true )
        .set { ch_genbank_acc_chunks }

    //// fetch Genbank sequences as .fasta
    FETCH_GENBANK (
        ch_genbank_acc_chunks,
        ch_entrez_key
    )

    //// add NCBI taxid to header of Genbank sequences 
    ADD_TAXID_GENBANK (
        FETCH_GENBANK.out.fetched_seqs
    )

    //// reformat sequence names to contain taoxnomic lineage
    RENAME_GENBANK (
        ADD_TAXID_GENBANK.out.fasta,
        GET_NCBI_TAXONOMY.out.rankedlineage_noname
    )

    //// populate ch_genbank_fasta channel
    RENAME_GENBANK.out.fasta
        .set { ch_genbank_fasta }

    //// count number of sequences downloaded from Genbank
    ch_count_genbank = ch_genbank_fasta.countFasta()
    
    /*
    Getting BOLD sequences and matching them to the NCBI database
    */

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

        //// chunk BOLD sequences into smaller .fasta files for processing
        MERGE_BOLD.out.fasta
            .splitText( by: params.chunk_size, file: true )
            .set { ch_bold_fasta }

        //// count number of sequences extracted from BOLD
        ch_count_bold = ch_bold_fasta.countFasta()

    } else {
        //// extract sequences by querying BOLD API
        // FETCH_BOLD (

        // )
        //// TODO: Need to make sure this outputs same tibble format as EXTRACT_BOLD
    }

    //// count number of internal sequences 
    ch_count_mito = ch_mito_fasta.countFasta()    

    //// count number of genome assembly-derived sequences 
    ch_count_genome = ch_genome_fasta.countFasta()

    // //// add internal sequences


    //// count number of internal sequences 
    ch_count_internal = ch_internal_fasta.countFasta()

    //// count number of external input sequences
    ch_count_external = ch_genbank_fasta
        .concat ( ch_bold_fasta )
        .concat ( ch_mito_fasta )
        .concat ( ch_genome_fasta )
        .countFasta() 

    //// concat output channels into a single channel for PHMM alignment
    ch_genbank_fasta
        .concat ( ch_bold_fasta )
        .concat ( ch_mito_fasta )
        .concat ( ch_genome_fasta )
        .concat ( ch_internal_fasta )
        .set { ch_input_seqs }  

    //// count number of total input sequences
    ch_count_input = ch_input_seqs.countFasta()      

    //// combine and save all unfiltered input sequences a single .fasta 
    if ( params.save_input || params.save_intermediate ) {
        ch_input_seqs
            .collectFile ( 
                name: "input_sequences.fasta",
                storeDir: "./output/results"
            )
    }

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
        ch_input_seqs,
        ch_phmm
    )

    //// combine and save intermediate file 
    if ( params.save_intermediate ) {
        FILTER_PHMM.out.fasta
            .collectFile ( 
                name: "filter_phmm.fasta",
                storeDir: "./output/results"
            )
    }

    //// count number of sequences passing PHMM filter
    ch_count_filter_phmm = FILTER_PHMM.out.fasta.countFasta() 

    //// optional: filter for stop codons
    if ( params.coding && params.genetic_code ){

        FILTER_STOP (
            FILTER_PHMM.out.seqs
        )

        //// combine and save intermediate file 
        if ( params.save_intermediate ) {
            FILTER_STOP.out.fasta
                .collectFile ( 
                    name: "filter_stop.fasta",
                    storeDir: "./output/results"
                )
        }

        //// count number of sequences passing stop codon filter
        ch_count_filter_stop = FILTER_STOP.out.fasta.countFasta() 

        ch_filter_output = FILTER_STOP.out.seqs.collect()

    } else {

        ch_filter_output = FILTER_PHMM.out.seqs.collect()

        ch_count_filter_stop = Channel.empty() // make empty channel as this filter isn't applied
    
    }

    // //// branch channels based on seq_source
    // ch_filter_output
    //     .branch {
    //         bold: it[1] == "bold"
    //         other: true
    //     }
    //     .set { ch_filter_branch }

    /// NOTE: Not currently needed as taxids are used
    // //// resolve taxonomic synonyms
    // RESOLVE_SYNONYMS ( 
    //     ch_filter_branch.other,
    //     GET_NCBI_TAXONOMY.out.db_path
    // )

    // //// collect chunks into a list
    // RESOLVE_SYNONYMS.out.seqs
    //     .concat ( ch_filter_branch.bold )
    //     .map { seqs, seq_source -> seqs } // remove seq_source
    //     .collect()
    //     .set { ch_chunks }

    //// combine sequences into one .fasta file
    COMBINE_CHUNKS ( 
        ch_filter_output
    )
 
    //// remove exact duplicate sequences if they exist
    REMOVE_EXACT_DUPLICATES (
        COMBINE_CHUNKS.out.fasta
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
    ch_count_remove_exact = REMOVE_EXACT_DUPLICATES.out.fasta.countFasta()

    //// remove contaminating sequences
    REMOVE_CONTAM (
        REMOVE_EXACT_DUPLICATES.out.fasta,
        GET_NCBI_TAXONOMY.out.rankedlineage
    )

    //// combine and save intermediate file 
    if ( params.save_intermediate ) {
        REMOVE_CONTAM.out.fasta
            .collectFile ( 
                name: "remove_contam.fasta",
                storeDir: "./output/results"
            )
    }

    //// count number of sequences passing taxonomic decontamination
    ch_count_remove_contam = REMOVE_CONTAM.out.fasta.countFasta()

    //// prune large groups
    PRUNE_GROUPS (
        REMOVE_CONTAM.out.fasta
    )

    //// count number of sequences passing group pruning
    ch_count_prune_groups = PRUNE_GROUPS.out.fasta.countFasta()

    // //// reformat names using taxonomic hierarchy
    // REFORMAT_NAMES (
    //     PRUNE_GROUPS.out.seqs,
    //     GET_NCBI_TAXONOMY.out.rankedlineage
    // )

    //// summarise number of taxa in database
    SUMMARISE_TAXA (
        PRUNE_GROUPS.out.seqs
    )
/* 
    //// summarise the counts of sequences at each stage of the pipeline
    SUMMARISE_COUNTS (
        ch_count_genbank,
        ch_count_bold,
        ch_count_mito,
        ch_count_genome,
        ch_count_internal,
        ch_count_external,
        ch_count_input,
        ch_count_filter_phmm,
        ch_count_filter_stop,
        ch_count_remove_exact,
        ch_count_remove_contam,
        ch_count_prune_groups
    )
 */

    //// sequence counts for debugging
    ch_count_genbank        .view { "ch_count_genbank: $it" }
    ch_count_bold           .view { "ch_count_bold: $it" }
    ch_count_mito           .view { "ch_count_mito: $it" }
    ch_count_genome         .view { "ch_count_genome: $it" }
    ch_count_internal       .view { "ch_count_internal: $it" }
    ch_count_external       .view { "ch_count_external: $it" }
    ch_count_input          .view { "ch_count_input: $it" }
    ch_count_filter_phmm    .view { "ch_count_filter_phmm: $it" }
    ch_count_filter_stop    .view { "ch_count_filter_stop: $it" }
    ch_count_remove_exact   .view { "ch_count_remove_exact: $it" }
    ch_count_remove_contam  .view { "ch_count_remove_contam: $it" }
    ch_count_prune_groups   .view { "ch_count_prune_groups: $it" }

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
    // curated_fasta = PRUNE_GROUPS.out.fasta
    // taxa_summary = TAXA_SUMMARY.out.csv
    // idtaxa_model = TRAIN_IDTAXA.out.model

}