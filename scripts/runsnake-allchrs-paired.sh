#!/bin/bash
NJOBS=50
#SPECIES="Bos_taurus"
#CHRS="1 2"
#SPECIES="Oryctolagus_cuniculus"
#SPECIES="Monodelphis_domestica"
#SPECIES="Equus_caballus"
#SPECIES="Canis_lupus_familiaris"
#CHRS="1 3"
#SPECIES="Capra_hircus"
#CHRS="18 19"
#SPECIES="Sus_scrofa"
#SPECIES="Ovis_aries"
#SPECIES="Bos_taurus"
#SPECIES="Macaca_mulatta"
SPECIES="Mus_musculus"
CHRS="3"
#CHRS="S 2 3"
PROJ="cmxci.sing.$SPECIES"
DATE=`date +"%Y%m%d.%H%M"`

echo 'Running on all chrs...' 
for chr in $CHRS; do 
	echo "Running on all samples against chromosome $chr"
	echo snakemake --keep-going --rerun-incomplete  --config species=$SPECIES chr=$chr --jobs $NJOBS --cluster "qsub -N $PROJ -pe threads {threads} -wd ~/work/cmxci.single -l m_mem_free={resources.mem_mb}M" 

	time snakemake --keep-going --rerun-incomplete  --config species=$SPECIES chr=$chr --jobs $NJOBS --cluster "qsub -N $PROJ -pe threads {threads} -wd ~/work/cmxci.single -l m_mem_free={resources.mem_mb}M"  > smlog.$DATE.txt 2>&1
done
