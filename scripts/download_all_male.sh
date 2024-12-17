#!/bin/bash
#SPECIES="Canis_lupus_familiaris Capra_hircus Equus_caballus Monodelphis_domestica Oryctolagus_cuniculus"
#SPECIES="Capra_hircus Equus_caballus Ovis_aries Rattus_norvegicus"
#SPECIES="Ovis_aries Rattus_norvegicus"
#SPECIES="Bos_taurus"
SPECIES="Sus_scrofa"
#SPECIES="Mus_musculus"
#SPECIES="Macaca_mulatta"
DATE=`date +"%Y%m%d.%H%M"`

for SPEC in $SPECIES; do 
	echo $SPEC
	mkdir -p $SPEC
        cat ${SPEC}_male_SRR.txt | xargs -I {} ~/git/cshlwork/scripts/sra_download.py -d -o $SPEC/ -r {} 	
done
