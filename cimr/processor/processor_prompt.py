#!/usr/bin/env python3

import logging
from .util import Infiler


def process_cli(args):
    outdir = args.outdir 
    outdir.mkdir(exist_ok=True)
    from .eqtl import checkgene
    from .tad import annotatetad
    """input files are checked for both type-dependant conditions and common formats"""
    if args.eqtl is not None:
        infile = Infiler('eqtl', args.eqtl)
        infile.readfile()
        checkgene()
    elif args.gwas is not None:
        infile = Infiler('gwas', args.gwas)
        infile.readfile()
    elif args.tad is not None:
        annotatetad(args.tad)
    else:
        logging.info(f'nothing to do')

