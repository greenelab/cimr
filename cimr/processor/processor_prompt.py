#!/usr/bin/env python3

import logging

def process_cli(args):
    outdir = args.outdir 
    outdir.mkdir(exist_ok=True)
    from cimr.processor.util import readfile
    from cimr.processor.eqtl import checkgene
    from cimr.processor.tad import annotatetad
    """input files are checked for both type-dependant conditions and common formats"""
    if args.eqtl is not None:
        readfile(args.eqtl)
    elif args.gwas is not None:
        readfile(args.gwas)
    elif args.tad is not None:
        annotatetad(args.tad)
    else:
        logging.info(f'nothing to do')

