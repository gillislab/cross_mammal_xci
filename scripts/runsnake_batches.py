#!/usr/bin/env python

import argparse
import datetime as dt
import logging
import os   
import sys
import traceback

import subprocess as sp

gitpath = os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from cshlwork.utils import *

#
#  MODE="single"
#  SPECIES="Ovis_aries"
#  CHR="3"
#  SAMPLE='["SRR11836547","SRR11836548"]'
#
def get_pathinfo(filepath):
    abspath = os.path.abspath(filepath)
    dirpath = os.path.dirname(abspath)
    file_name = os.path.basename(abspath)
    (base, ext) = os.path.splitext(file_name)
    ext = ext[1:] # remove dot
    return (dir, base, ext)
    

def run_snake_batch(mode, samples, species, chromosome, numjobs, latency ):
    STATUS="/grid/gillis/home/hover/git/elzar-example/snakemake/status-sge.py"
    proj=f"cmxci.{mode}.{species}"
    datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    sampstr = str(samples).replace(' ','')
    

    cmd = [ 
            'snakemake',
            '--keep-going',
            '--rerun-incomplete',
            '--latency-wait', str(latency) ,
            '--config', f'species={species}',f'chr={chromosome}',f'sample={sampstr}', 
            '--jobs', str(numjobs),
            '--cluster-status', STATUS,
            '--cluster-cancel','qdel',
            '--cluster', f'qsub -N {proj} -pe threads {{threads}} -wd /grid/gillis/home/hover/work/cmxci.{mode} -l m_mem_free={{resources.mem_mb}}M '
          ]
    try:
        cmdstr = ' '.join(cmd)
        logging.info(f'    Running command: {cmdstr}')
        cp = run_command(cmd)
        logging.info(f'got rc={cp.returncode}')
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
                        required=True,
                        default=None, 
                        help='species to run on. ')

    parser.add_argument('-c', '--chromosome', 
                        metavar='chromosome', 
                        type=str,
                        default='1', 
                        help='chromsome ')

    parser.add_argument('-j','--numjobs', 
                        metavar='numjobs',
                        required=False,
                        default=2,
                        type=int, 
                        help='Max simultaneous jobs.')

    parser.add_argument('-b','--batchsize', 
                        metavar='batchsize',
                        required=False,
                        default=2,
                        type=int, 
                        help='Number of samples per batch.')

    parser.add_argument('-L','--latency', 
                        metavar='latency',
                        required=False,
                        default=15,
                        type=int, 
                        help='Seconds to wait for output files.')

    parser.add_argument('infiles', 
                        metavar='infiles', 
                        nargs='+',
                        type=str, 
                        help='SRA files to process. ') 
   
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    logging.debug(args)


    samplelist = []
    for filename in args.infiles:
        (dir, base, ext) = get_pathinfo(filename)
        samplelist.append(base)
    
    batchlist = [samplelist[i:i + args.batchsize] for i in range(0, len(samplelist), args.batchsize)]
    logging.debug(f'batchlist= {batchlist}')
    
    batchno = 0
    for batch in batchlist:
        batchno += 1
        logging.debug(f'executing w/ mode={args.mode} species={args.species} chr={args.chromosome} numjobs={args.numjobs}')
        logging.info(f'executing batch [{batchno}/{len(batchlist)}] w/ samples={batch}')
        run_snake_batch(args.mode, batch, args.species, args.chromosome, args.numjobs, args.latency )
        logging.info(f'finished invocation of snakemake...') 
    
      