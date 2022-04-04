#!/bin/bash
NJOBS=20
SPECIES="Sus_scrofa"
#SPECIES="Ovis_aries"
#SPECIES="Bos_taurus"
#SPECIES="Macaca_mulatta"
PROJ="cmxci.paired.$SPECIES"
DATE=`date +"%Y%m%d.%H%M"`

echo 'Running on all samples...'
echo snakemake --keep-going --config species=$SPECIES --jobs $NJOBS --cluster "qsub -N $PROJ -pe threads {threads} -wd ~/work/cmxci.paired -l m_mem_free={resources.mem_mb}M" 

time snakemake --keep-going --config species=$SPECIES --jobs $NJOBS --cluster "qsub -N $PROJ -pe threads {threads} -wd ~/work/cmxci.paired -l m_mem_free={resources.mem_mb}M" > smlog.$DATE.txt 2>&1 &

