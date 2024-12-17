#!/bin/bash
NJOBS=30
#SPECIES="Bos_taurus"
#SPECIES="Oryctolagus_cuniculus"
#SPECIES="Monodelphis_domestica"
#SPECIES="Sus_scrofa"
#SPECIES="Canis_lupus_familiaris"
#SPECIES="Capra_hircus"
#SPECIES="Equus_caballus"
MODE="single"
#SPECIES="Ovis_aries"
SPECIES="Rattus_norvegicus"
#SPECIES="Bos_taurus"
#SPECIES="Macaca_mulatta"
PROJ="cmxci.$MODE.$SPECIES"
DATE=`date +"%Y%m%d.%H%M"`
CHR="3"

echo 'Running snakemake'
#echo snakemake --keep-going --config species=$SPECIES --jobs $NJOBS --cluster "qsub -N $PROJ -pe threads {threads} -wd ~/work/cmxci.$MODE -l m_mem_free={resources.mem_mb}M" 
#time snakemake --keep-going --config species=$SPECIES --jobs $NJOBS --cluster "qsub -N $PROJ -pe threads {threads} -wd ~/work/cmxci.$MODE -l m_mem_free={resources.mem_mb}M" > smlog.$DATE.txt 2>&1 &

echo snakemake  --config species=$SPECIES chr=$CHR --jobs $NJOBS --executor cluster-generic --cluster-generic-submit-cmd \"qsub -N $PROJ -pe threads {threads} -wd ~/work/cmxci.$MODE -l m_mem_free={resources.mem_mb}M\" 
#time snakemake  --config species=$SPECIES chr=$CHR --jobs $NJOBS --executor cluster-generic --cluster-generic-submit-cmd "qsub -N $PROJ -pe threads {threads} -wd ~/work/cmxci.$MODE -l m_mem_free={resources.mem_mb}M" 

