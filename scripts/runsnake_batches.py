#/usr/bin/env python

import argparse
import datetime as dt
import logging
import os   
import traceback

import subprocess as sp

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from cshlwork.utils import *

#
#  MODE="single"
#  SPECIES="Ovis_aries"
#  CHR="3"
#  SAMPLE='["SRR11836547","SRR11836548"]'
#
#






MODE="paired"
PROJ=f"cmxci.{MODE}"
DATESTR = dt.datetime.now().strftime("%Y%m%d%H%M")
STATUS="/grid/gillis/home/hover/git/elzar-example/snakemake/status-sge.py"

logging.getLogger().setLevel(logging.DEBUG)
logging.debug(f"PATH={os.environ['PATH']}")

for SPEC in SPECIES.keys():
    logging.info(f'Running on {SPEC}')
    if FIRSTONLY:
        logging.info(f'only running first chromosome')
        CHRS = []
        CHRS.append(SPECIES[SPEC][0])
    else:
        CHRS = SPECIES[SPEC]

    for CHR in CHRS:
        logging.info(f'    Running {SPEC} -> {CHR}')
        cmd = [ 
                'snakemake',
                '--keep-going',
                '--rerun-incomplete',
                '--latency-wait', '10',
                '--config',f'species={SPEC}',f'chr={CHR}',
                '--jobs', str(NJOBS),
                '--cluster-status', STATUS,
                '--cluster-cancel','qdel',
                '--cluster', f'qsub -N {PROJ} -pe threads {{threads}} -wd /grid/gillis/home/hover/work/cmxci.{MODE} -l m_mem_free={{resources.mem_mb}}M'
              ]
        try:
            cmdstr = ' '.join(cmd)
            logging.info(f'    Running command: {cmdstr}')
            run_command(cmd)
        except NonZeroReturnException as nzre:
            logging.error(traceback.format_exc(None))





if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARN)
    
    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')

    parser.add_argument('-m','--mode', 
                        metavar='mode', 
                        type=str,
                        default='paired', 
                        help='mode to run in [single|paired] ')
    
    parser.add_argument('-s', '--species', 
                        metavar='species', 
                        type=str,
                        default=None, 
                        help='species to run on. ')
    
    parser.add_argument('-j','--numjobs', 
                        metavar='numjobs',
                        required=False,
                        default=5,
                        type=int, 
                        help='Max simultaneous jobs.')

    parser.add_argument('-b','--batchsize', 
                        metavar='batchsize',
                        required=False,
                        default=5,
                        type=int, 
                        help='Number of samples per batch.')

   
#    parser.add_argument('srafiles', 
#                        metavar='srafiles', 
#                        nargs='+',
#                        type=argparse.FileType('r'), 
#                        help='SRA files to process. ')
#
#        for f in args.srafiles:
#               for line in f:
#                      process...
#   

    parser.add_argument('srafiles', 
                        metavar='srafiles', 
                        nargs='+',
                        type=str, 
                        help='SRA files to process. ') 
   
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    logging.debug(args)











NJOBS=2
#SPECIES="Bos_taurus"
#CHR="1"

#SPECIES="Canis_lupus_familiaris"

#SPECIES="Capra_hircus"

#SPECIES="Equus_caballus"


#SPECIES="Macaca_mulatta"
#SPECIES="Monodelphis_domestica"
#SPECIES="Oryctolagus_cuniculus"
#SPECIES="Sus_scrofa"



PROJ="cmxci.$MODE.$SPECIES"
DATE=`date +"%Y%m%d.%H%M"`

echo 'Running on all samples...'
echo snakemake --keep-going --config species=$SPECIES sample=$SAMPLE chr=$CHR --jobs $NJOBS --cluster "qsub -N $PROJ -pe threads {threads} -wd ~/work/cmxci.$MODE -l m_mem_free={resources.mem_mb}M" 

time snakemake --keep-going --config species=$SPECIES sample=$SAMPLE chr=$CHR --jobs $NJOBS --cluster "qsub -N $PROJ -pe threads {threads} -wd ~/work/cmxci.$MODE -l m_mem_free={resources.mem_mb}M" > smlog.$DATE.txt 2>&1