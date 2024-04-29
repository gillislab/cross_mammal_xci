#!/usr/bin/env python
#
#  Check SRA files to determine if they used paired-end or single-end sequencing. 
#
# See documentation for the SRA Toolkit:
# http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump
#
#
# 

import argparse
import logging
import os   
import traceback

import subprocess as sp

def handle_file(fn, single=None, paired=None):
    basename = os.path.basename(fn)
    dirname = os.path.dirname(fn)
    
    if single is not None:
        single = os.path.abspath( os.path.expanduser(single) )
    if paired is not None:
        paired = os.path.abspath( os.path.expanduser(paired) )
       
    try:
        b = is_paired(fn)
        logging.debug(f'{fn} -> {b}')
        if b:
            print(f'{fn} is paired.')
        else:
            print(f'{fn} is single.')
        
        if single is not None:
            if not b:
                os.makedirs(f'{single}', exist_ok = True)
                logging.info(f'{fn} -> {single}/{basename}')
                os.rename( fn , f'{single}/{basename}' )
        
        if paired is not None:
            if b:
                os.makedirs(f'{paired}', exist_ok = True)
                logging.info(f'{fn} -> {paired}/{basename}')
                os.rename( fn , f'{paired}/{basename}' )                    
               
    except Exception as e:
        logging.error(f'problem with {fn}')
        logging.error(traceback.format_exc(None))


def is_paired(filename):
    filename = os.path.abspath(filename)
    try:
        contents = sp.check_output(["fastq-dump","-X","1", "-Z", "--split-spot", filename],
                                   stderr = sp.DEVNULL
                                   )
        logging.debug(f'contents = {contents}')
    except sp.CalledProcessError as e:
        raise Exception(f"Error running fastq-dump on {filename}");

    if(contents.count(b"\n") == 4):
        return False

    elif(contents.count(b"\n") == 8):
        return True
    else:
        raise Exception("Unexpected output from fastq-dump on ", filename)
    


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

    parser.add_argument('-s','--single', 
                        metavar='single', 
                        type=str,
                        default=None, 
                        help='directory to move single-ended to')
    
    parser.add_argument('-p', '--paired', 
                        metavar='paired', 
                        type=str,
                        default=None, 
                        help='directory to move paired-ended to')
        
   
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

    total_files = len(args.srafiles)
    handled_files = 0    
    for fn in args.srafiles:
        logging.info(f'handling {handled_files + 1} of {total_files}...')
        handle_file(fn, args.single, args.paired)
        handled_files += 1
       

    
    