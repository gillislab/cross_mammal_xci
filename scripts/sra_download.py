#!/usr/bin/env python
import argparse
import os
import sys
import pandas as pd
import logging

from configparser import ConfigParser

gitpath = os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from cshlwork.utils import *
from sra.utils import *

if __name__ == "__main__":
    FORMAT = '%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
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

    parser.add_argument('-c', '--config',
                            action="store",
                            dest='conffile',
                            default=None,
                            help='Config file path [~/git/cshlwork/etc/sra.conf]')

    parser.add_argument('-f', '--fasterq',
                        metavar='fastqrun',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download args with fasterq-dump. e.g. SRR14584407')

    parser.add_argument('-p', '--projects',
                        metavar='proj_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download project files, e.g. SRP151064')

    parser.add_argument('-e', '--experiments',
                        metavar='exp_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download experiment files, e.g. SRX2798918')


    parser.add_argument('-r', '--runs',
                        metavar='run_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download run SRA file, e.g. SRR14584407')

    parser.add_argument('-s', '--samples',
                        metavar='samp_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download project files, e.g. SRP151064')

    parser.add_argument('-o', '--outdir',
                        metavar='outdir',
                        type=str,
                        default=None,
                        help='Outdir [ <cwd> ] ')

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    if args.conffile is not None:
        cp = ConfigParser()
        cp.read(os.path.expanduser(args.conffile)) 
    else:
        cp = get_default_config()
        
    cs = get_configstr(cp)
    logging.debug(f"got config: {cs}")

    logging.debug(f"args: {args}")
    
    if args.runs is not None:
        for run_id in args.runs:
            logging.debug(f'Downloading SRA for run {run_id}')
            download_run_sra(run_id, args.outdir)
        logging.info('Done')


    if args.samples is not None:
        for samp_id in args.samples:
            logging.debug(f'Querying SRA for sample {samp_id}')
            #download_run_sra(run_id, args.outdir)
        logging.info('Done')

    if args.projects is not None:
        for proj_id in args.projects:
            logging.debug(f'Querying SRA for project {proj_id}')
            df = query_project_metadata(proj_id)
            print(df)
            #download_run_sra(run_id, args.outdir)
        logging.info('Done')        


