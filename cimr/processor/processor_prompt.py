#!/usr/bin/env python3


import logging
from .utils import Infiler
from .utils import Integrator
from .query import Querier
from .tad import annotate_tad
    

def processor_cli(args):
    
    data_types = {'gwas', 'eqtl'}
    annotations = {'tad'}
    outdir = args.outdir 
    outdir.mkdir(exist_ok=True)
    logging.info(f' directory {str(outdir)} will be used for cimr jobs.')
    outfile = str(outdir) + '/' + str(args.out) + '_'
    logging.info(f' file prefix {str(args.out)} will be used for output files')

    data_type = args.data_type

    if data_type in data_types:
        if args.file_name is not None:
            outfile = outfile + data_type + '.txt'
            infile = Infiler(data_type, args.file_name, args.genome_build)
            infile.read_file()
            genes = list(infile.list_genes())
            if data_type == 'eqtl':
                Querier(genes)
            infile.write_file(outfile)

    elif data_type in annotations:
        annotate_tad(args.file_name)

    else:
        logging.error(f' data_type or file_name is not recognized. nothing to do.')
    
    if args.process:
        logging.info(f' in order to contribute to the cimr database, use --integrate option.')
    elif args.integrate:
        if data_type in data_types:
            integrating = Integrator(
                data_type, 
                args.file_name, 
                can_be_public=args.can_be_public, 
                genome_build=args.genome_build
            )
            integrating.make_local_db(args.temp_dir)
    else:
        logging.error(f' did the command include either --process or --integrate functions?')
        

