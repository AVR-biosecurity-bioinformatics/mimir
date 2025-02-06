#!/usr/bin/env bash
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
        TYPE="mitochondrial"
        PHMM_URL="https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00115?annotation=hmm"
        SEED_URL="https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00115/?annotation=alignment:seed&download"
        ;;
    
    # "28S")
    #     GENBANK="28S[TI] OR LSU[TI] OR large ribosomal subunit[TI]"
    #     BOLD="28S"
    #     CODING="false"
    #     TYPE="nuclear"
    #     ;;
    
    #"ITS")
    #"matK")
    #"rbcL")

    *)
        echo "Marker value '${3}' not valid"
        exit 1
        ;;

esac

### TODO: make these downloads more robust

# download PHMM from InterPro and unzip
curl $PHMM_URL | gunzip > full.hmm

# download seed alignment for PHMM from InterPro and unzip
curl $SEED_URL | gunzip > seed.stockholm