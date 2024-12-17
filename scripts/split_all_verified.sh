#!/bin/bash
#SPECIES="Canis_lupus_familiaris Capra_hircus Equus_caballus Monodelphis_domestica Oryctolagus_cuniculus"
#SPECIES="Bos_taurus Canis_lupus_familiaris Capra_hircus Equus_caballus Macaca_mulatta Ovis_aries Rattus_norvegicus Sus_scrofa"
#SPECIES="Mus_musculus"
#SPECIES="Macaca_mulatta"
#SPECIES="Rattus_norvegicus"
#SPECIES="Bos_taurus"
SPECIES="Sus_scrofa"
DATE=`date +"%Y%m%d.%H%M"`

for SPEC in $SPECIES; do 
	echo $SPEC
	mkdir -p $SPEC.paired $SPEC.single
        time ~/git/cross_mammal_xci/scripts/sra-paired.py -v -s $SPEC.single -p $SPEC.paired $SPEC/*.sra 	
done
