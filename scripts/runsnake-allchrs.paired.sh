#!/bin/bash
NJOBS=40
SPECIES="Bos_taurus"
#CHRS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 X"
#SPECIES="Oryctolagus_cuniculus"
#SPECIES="Monodelphis_domestica"
#SPECIES="Equus_caballus"
#SPECIES="Canis_lupus_familiaris"
CHRS="1 3"
#SPECIES="Capra_hircus"
#SPECIES="Sus_scrofa"
#SPECIES="Ovis_aries"
#SPECIES="Bos_taurus"
#SPECIES="Macaca_mulatta"
PROJ="cmxci.paired.$SPECIES"
DATE=`date +"%Y%m%d.%H%M"`

echo 'Running on all chrs...' 
for chr in $CHRS; do 
	echo "Running on all samples against chromosome $chr"
	echo snakemake --keep-going --rerun-incomplete  --config species=$SPECIES chr=$chr --jobs $NJOBS --cluster "qsub -N $PROJ -pe threads {threads} -wd ~/work/cmxci.paired -l m_mem_free={resources.mem_mb}M" 

	time snakemake --keep-going --rerun-incomplete  --config species=$SPECIES chr=$chr --jobs $NJOBS --cluster "qsub -N $PROJ -pe threads {threads} -wd ~/work/cmxci.paired -l m_mem_free={resources.mem_mb}M" > smlog.$DATE.txt 2>&1
done
