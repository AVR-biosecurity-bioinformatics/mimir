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

if( !nextflow.version.matches('=23.04.5') ) {
    println " "
    println "*** ERROR ~ This pipeline currently requires Nextflow version 23.04.5 -- You are running version ${nextflow.version}. ***"
    error "*** You can use version 23.04.5 by appending 'NXF_VER=23.04.5' to the front of the 'nextflow run' command. ***"
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
validateParameters( parameters_schema: 'nextflow_schema.json' )


// Print summary of supplied parameters (that differ from defaults)
log.info paramsSummaryLog(workflow)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//// define main workflow
workflow MIMIR {

    println "Hello"

}

///// run implicit workflow
workflow {
    MIMIR()
}