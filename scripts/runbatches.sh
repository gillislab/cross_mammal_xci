#!/bin/bash
# time snakemake --keep-going --batch all=1/5 --config species=Bos_taurus --jobs 20 --cluster "qsub -N cmxci1 -pe threads {threads} -wd ~/work/cmxci1 -l m_mem_free={resources.mem_mb}M" > smlog.3.txt 2>&1

#BATCHES="1 2 3 4 5 6 7 8 9 10"
BATCHES="1"
#BATCHES="5 8 9"
TOTAL=4
NJOBS=4
SPECIES="Sus_scrofa"
#SPECIES="Ovis_aries"
#SPECIES="Bos_taurus"
#SPECIES="Macaca_mulatta"
PROJ="cmxci.paired"
DATE=`date +"%Y%m%d.%H%M"`
WORKFLOW="/grid/gillis/home/hover/git/cross_mammal_xci/workflow.paired"

for n in $BATCHES; do
	mkdir -p ./batch$n
	cd ./batch$n
	ln -s $WORKFLOW workflow
	cd ..
done


for n in $BATCHES; do
	echo time snakemake --keep-going --rerun-incomplete -d ./batch$n --batch all=$n/$TOTAL --config species=$SPECIES --jobs $NJOBS --cluster \"qsub -N $PROJ -pe threads {threads} -wd ~/work/$PROJ -l m_mem_free={resources.mem_mb}M\" > cmd.$DATE.$n.txt 2>&1
	time snakemake --keep-going --rerun-incomplete -d ./batch$n --batch all=$n/$TOTAL --config species=$SPECIES --jobs $NJOBS --cluster "qsub -N $PROJ -pe threads {threads} -wd ~/work/$PROJ -l m_mem_free={resources.mem_mb}M" > smlog.$DATE.$n.txt 2>&1 & 
done

