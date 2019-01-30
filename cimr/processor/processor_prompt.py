#!/usr/bin/env python3

import logging
from .utils import Infiler
from .utils import Integrator
from .query import Querier
from .tad import annotate_tad
    

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
        if args.filename is not None:
            outfile = outfile + datatype + '.txt'
            infile = Infiler(datatype, args.filename, args.genome_build)
            infile.read_file()
            genes = infile.list_genes()
            if datatype == 'eqtl':
                Querier(genes)
            infile.write_file(outfile)

    elif datatype in annotations:
        annotate_tad(args.filename)

    else:
        logging.error(f' datatype or filename is not recognized. nothing to do.')
    
    if args.process:
        logging.info(f' in order to contribute to the cimr database, use --integrate option.')
    elif args.integrate:
        if datatype in datatypes:
            integrating = Integrator(
                datatype, 
                args.filename, 
                can_be_public=args.can_be_public, 
                genome_build=args.genome_build
            )
            integrating.make_local_db(args.tempdir)
    else:
        logging.error(f' did the command include either --process or --integrate functions?')
        

