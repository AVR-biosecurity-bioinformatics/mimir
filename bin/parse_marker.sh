#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = cpus
# $3 = marker

case ${3} in 
    
    "COI")
        GENBANK="COI[GENE] OR COX1[GENE] OR COXI[GENE]"
        BOLD="COI-5P"
        CODING="true"
        ;;
    
    "28S")
        GENBANK="28S[TI] OR LSU[TI] OR large ribosomal subunit[TI]"
        BOLD="28S"
        CODING="false"
        ;;
    
    #"ITS")
    #"matK")
    #"rbcL")

    *)
        echo "Marker value '${3}' not valid"
        exit 1
        ;;

esac
