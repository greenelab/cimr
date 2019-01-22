#!/usr/bin/env python3

import logging
from .utils import Infiler
from .utils import Integrator
from .eqtl import checkgene
from .tad import annotatetad
    

def processor_cli(args):
    
    datatypes = {'gwas', 'eqtl'}
    annotations = {'tad'}
    outdir = args.outdir 
    outdir.mkdir(exist_ok=True)
    logging.info(f' directory {str(outdir)} will be used for cimr jobs.')
    outfile = str(outdir) + '/' + str(args.out) + '_'
    logging.info(f' file prefix {str(args.out)} will be used for output files')

    datatype = args.datatype

    if datatype in datatypes:
        outfile = outfile + datatype + '.txt'
        infile = Infiler(datatype, args.filename, args.genome_build)
        infile.readfile()
        if datatype == 'eqtl':
            checkgene()
        infile.writefile(outfile)

    elif datatype in annotations:
        annotatetad(args.filename)

    else:
        logging.error(f' datatype or filename is not recognized. nothing to do.')
    
    if process:
        logging.info(f' in order to contribute to the cimr database, use --integrate option.')
    elif integrate:
        if datatype in datatypes:
            

