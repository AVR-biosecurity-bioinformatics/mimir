#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AVR-biosecurity-bioinformatics/mimir
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/AVR-biosecurity-bioinformatics/mimir
----------------------------------------------------------------------------------------
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include functions from nf-schema
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema' 

// 
// from https://github.com/nextflow-io/nextflow/issues/1129

// if( !nextflow.version.matches('=23.04.5') ) {
//     println " "
//     println "*** ERROR ~ This pipeline currently requires Nextflow version 23.04.5 -- You are running version ${nextflow.version}. ***"
//     error "*** You can use version 23.04.5 by appending 'NXF_VER=23.04.5' to the front of the 'nextflow run' command. ***"
// }

if( !nextflow.version.matches('=23.05.0-edge') ) {
    println " "
    println "*** ERROR ~ This pipeline currently requires Nextflow version 23.05.0-edge -- You are running version ${nextflow.version}. ***"
    error "*** You can use version 23.05.0-edge by appending 'NXF_VER=23.05.0-edge' to the front of the 'nextflow run' command. ***"
}

def startupMessage() {
    log.info pipelineHeader()
    log.info "~~~ mimir: DNA barcode reference database curation ~~~"
    log.info " "
}

def pipelineHeader(){
    return """                                                                          
                     ##               ##       
                                                    
          #####     ###    #####     ###    ###  ## 
          # # #     ##     # # #     ##     ###   
          # # #     #      # # #     #      ##    
         # # #     ##     # # #     ##     ##     
        ## # ##  ######  ## # ##  ######  #### 
    """.stripIndent()
}

startupMessage()

workflow.onComplete {
    if ( workflow.success ) {
      log.info "[$workflow.complete] >> Pipeline finished SUCCESSFULLY after $workflow.duration"
    } else {
      log.info "[$workflow.complete] >> Pipeline finished with ERRORS after $workflow.duration"
    }
    /*
    TODO: use other metadata (https://www.nextflow.io/docs/latest/metadata.html) to display at the end of a run.
    */

}


// Print help message, supply typical command line usage for the pipeline
if (params.help) {
//    log.info startupMessage()
   log.info paramsHelp("nextflow run AVR-biosecurity-bioinformatics/mimir") // TODO: add typical commands for pipeline
   exit 0
}

// Validate input parameters using schema
// validateParameters( parameters_schema: 'nextflow_schema.json' )


// Print summary of supplied parameters (that differ from defaults)
log.info paramsSummaryLog(workflow)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//// subworkflows
include { FILTER_SEQUENCES                                             } from './nextflow/subworkflows/filter_sequences'
include { FORMAT_DATABASE                                             } from './nextflow/subworkflows/format_database'
include { GET_BOLD                                             } from './nextflow/subworkflows/get_bold'
include { GET_GENBANK                                             } from './nextflow/subworkflows/get_genbank'
include { GET_INTERNAL                                             } from './nextflow/subworkflows/get_internal'
include { PARSE_INPUTS                                             } from './nextflow/subworkflows/parse_inputs'
include { SUMMARISE_DATABASE                                             } from './nextflow/subworkflows/summarise_database'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//// define main workflow
workflow MIMIR {

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


    /*
    Parse pipeline inputs and fetch taxonomic and marker-related information
    */

    PARSE_INPUTS (
        "dummy",
        ch_entrez_key,
        ch_targets,
        ch_marker
    )

    /*
    Get sequences from GenBank 
    */

    if ( params.use_genbank ){

        GET_GENBANK (
            PARSE_INPUTS.out.ch_taxon_idrank,
            PARSE_INPUTS.out.ch_genbank_query,
            ch_entrez_key,
            PARSE_INPUTS.out.ch_rankedlineage_noname
        )

        ch_genbank_fasta = GET_GENBANK.out.ch_genbank_fasta
        
    } else {
        
        ch_genbank_fasta = Channel.empty()
    
    }

    //// collect GenBank sequences into a single .fasta for source tracking
    ch_genbank_fasta
        .collectFile ( name: 'genbank.fasta', newLine: true, cache: true )
        .set { ch_source_genbank }

    //// count sequences from GenBank
    ch_count_genbank = ch_source_genbank.countFasta().combine(["genbank"])


    /*
    Get sequences from BOLD
    */

    if ( ( params.bold_db_path || params.bold_db_url ) && params.use_bold ) {

        GET_BOLD (
            ch_bold_db_path,
            ch_bold_db_url,
            PARSE_INPUTS.out.ch_bold_targets,
            PARSE_INPUTS.out.ch_bold_query,
            PARSE_INPUTS.out.ch_lineageparents,
            PARSE_INPUTS.out.ch_synonyms,

        )

        ch_bold_fasta = GET_BOLD.out.ch_bold_fasta

    } else {
    
        ch_bold_fasta = Channel.empty()
    
    }

    //// collect BOLD sequences into a single .fasta for source tracking
    ch_bold_fasta
        .collectFile ( name: 'bold.fasta', newLine: true, cache: true )
        .set { ch_source_bold }

    //// count number of BOLD sequences
    ch_count_bold = ch_bold_fasta.countFasta().combine(["bold"])


    /*
    Import and format internal (user-supplied) sequences
    */

    if ( ch_internal_seqs ){

        GET_INTERNAL (
            ch_internal_seqs
        )

        ch_internal_fasta = GET_INTERNAL.out.ch_internal_fasta
        ch_internal_names = GET_INTERNAL.out.ch_internal_names

    } else {

        ch_internal_fasta = Channel.empty()
        ch_internal_names = Channel.empty()

    }

    //// collect internal sequences into a single .fasta for source tracking
    GET_INTERNAL.out.ch_internal_fasta
        .collectFile ( name: 'internal.fasta', newLine: true, cache: true )
        .set { ch_source_internal }

    //// count number of internal sequences 
    ch_count_internal = ch_internal_fasta.countFasta().combine(["internal"])


    /*
    Get mitochondrial genomes from NCBI
    */

    //// PLACEHOLDER for dedicated mito subworkflow
    ch_mito_fasta               = Channel.empty()

    //// collect mitochondrial genome sequences into a single .fasta for source tracking
    ch_mito_fasta
        .collectFile ( name: 'mito.fasta', newLine: true, cache: true )
        .set { ch_source_mito }

    //// count number of mitochondrial sequences 
    ch_count_mito = ch_mito_fasta.countFasta().combine(["mito"])


    /*
    Get whole genomes from NCBI
    */

    //// PLACEHOLDER for dedicated genome subworkflow
    ch_genome_fasta             = Channel.empty()

    //// count number of genome assembly-derived sequences 
    ch_count_genome = ch_genome_fasta.countFasta().combine(["genome"])

    //// collect whole genome sequences into a single .fasta for source tracking
    ch_genome_fasta
        .collectFile ( name: 'genome.fasta', newLine: true, cache: true )
        .set { ch_source_genome }


    /*
    Combined input channels 
    */

    //// count number of external input sequences
    ch_count_external = ch_genbank_fasta
        .mix ( ch_bold_fasta )
        .mix ( ch_mito_fasta )
        .mix ( ch_genome_fasta )
        .countFasta()
        .combine(["external"])

    //// concat output channels into a single channel for PHMM alignment
    Channel.empty()
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


    /*
    Filter sequences
    */

    FILTER_SEQUENCES (
        ch_input_seqs,
        ch_internal_names,
        PARSE_INPUTS.out.ch_marker_coding,
        PARSE_INPUTS.out.ch_marker_type,
        PARSE_INPUTS.out.ch_gencodes,
        PARSE_INPUTS.out.ch_full_phmm,
        PARSE_INPUTS.out.ch_trimmed_phmm,
        PARSE_INPUTS.out.ch_frame_info
    )


    /*
    Format final database
    */

    FORMAT_DATABASE (
        FILTER_SEQUENCES.out.ch_final_seqs
    )


    /*
    Summarise database
    */

    //// collect sequence origins into a list of .fasta files for conversion into .csv
    Channel.empty()
        .concat ( ch_source_genbank )
        .concat ( ch_source_bold )
        .concat ( ch_source_mito )
        .concat ( ch_source_genome )
        .concat ( ch_source_internal )
        .collect ()
        .set { ch_sources }

    SUMMARISE_DATABASE (
        FORMAT_DATABASE.out.ch_final_database,
        ch_sources,
        FILTER_SEQUENCES.out.ch_fates,
        ch_key_species_list
    )


    /*
    Debugging code
    */

    //// sequence counts for debugging
    if ( params.debug_mode ){
        ch_count_genbank                                        .view{ "${it[1]}: ${it[0]}" }
        ch_count_bold                                           .view{ "${it[1]}: ${it[0]}" }
        ch_count_mito                                           .view{ "${it[1]}: ${it[0]}" }
        ch_count_genome                                         .view{ "${it[1]}: ${it[0]}" }
        ch_count_internal                                       .view{ "${it[1]}: ${it[0]}" }
        ch_count_external                                       .view{ "${it[1]}: ${it[0]}" }
        ch_count_input                                          .view{ "${it[1]}: ${it[0]}" }
        FILTER_SEQUENCES.out.ch_count_filter_unclassified       .view{ "${it[1]}: ${it[0]}" }
        FILTER_SEQUENCES.out.ch_count_filter_phmm_full          .view{ "${it[1]}: ${it[0]}" }
        FILTER_SEQUENCES.out.ch_count_filter_phmm_trimmed       .view{ "${it[1]}: ${it[0]}" }
        FILTER_SEQUENCES.out.ch_count_filter_duplicates         .view{ "${it[1]}: ${it[0]}" }
        FILTER_SEQUENCES.out.ch_count_filter_ambiguous          .view{ "${it[1]}: ${it[0]}" }
        FILTER_SEQUENCES.out.ch_count_filter_tax_outliers       .view{ "${it[1]}: ${it[0]}" }
        FILTER_SEQUENCES.out.ch_count_filter_redundant          .view{ "${it[1]}: ${it[0]}" }
        FILTER_SEQUENCES.out.ch_count_filter_seq_outliers       .view{ "${it[1]}: ${it[0]}" }
        FILTER_SEQUENCES.out.ch_count_select_final_sequences    .view{ "${it[1]}: ${it[0]}" }
    }


}

///// run implicit workflow
workflow {
    MIMIR()
}