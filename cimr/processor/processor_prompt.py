#!/usr/bin/env python3

import logging

def processor_cli(args):
    from cimr.processor.util import readfile
    from cimr.processor.util import writefile
    from cimr.processor.eqtl import checkgene
    from cimr.processor.tad import annotatetad
    
    outdir = args.outdir 
    outdir.mkdir(exist_ok=True)
    logging.info(f' directory {str(outdir)} will be used for cimr jobs.')

    outfile = str(outdir) + '/' + str(args.out) + '_'
    logging.info(f' file prefix {str(args.out)} will be used for output files')
    
    """input files are checked for both type-dependant conditions and common formats"""
    if args.eqtlfile is not None:
        outfile = outfile + 'eqtl.txt'
        summary_data = readfile(args.eqtl)
        summary_data = checkgene(summary_data)
        writefile(summary_data, outfile)
    elif args.gwasfile is not None:
        outfile = outfile + 'gwas.txt'
        summary_data = readfile(args.gwasfile)
        writefile(summary_data, outfile)
    elif args.tadfile is not None:
        outfile = outfile + 'tad.txt'
        annotatetad(args.tad)
    else:
        logging.info(f' nothing to do')

