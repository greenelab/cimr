#!/usr/bin/env python3

import logging
from .utils import Infiler
from .eqtl import checkgene
from .tad import annotatetad
    

def processor_cli(args):
    outdir = args.outdir 
    outdir.mkdir(exist_ok=True)
    logging.info(f' directory {str(outdir)} will be used for cimr jobs.')
    outfile = str(outdir) + '/' + str(args.out) + '_'
    logging.info(f' file prefix {str(args.out)} will be used for output files')

    """input files are checked for both type-dependant conditions and common formats"""
    if args.eqtlfile is not None:
        outfile = outfile + 'eqtl.txt'
        infile = Infiler('eqtl', args.eqtlfile, args.genome_build)
        checkgene()
        infile.readfile()
        infile.writefile(outfile)
    elif args.gwasfile is not None:
        outfile = outfile + 'gwas.txt'
        infile = Infiler('gwas', args.gwasfile, args.genome_build)
        infile.readfile()
        infile.writefile(outfile)
    elif args.tadfile is not None:
        outfile = outfile + 'tad.txt'
        annotatetad(args.tadfile)
    else:
        logging.info(f' nothing to do.')

