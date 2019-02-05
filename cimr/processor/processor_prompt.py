#!/usr/bin/env python3


import logging
from .utils import Infiler
from .utils import Integrator
from .query import Querier
from .tad import annotate_tad
    

def processor_cli(args):
    
    data_types = {'gwas', 'eqtl'}
    annotations = {'gene', 'tad'}
    outdir = args.outdir 
    outdir.mkdir(exist_ok=True)
    logging.info(f' directory {str(outdir)} will be used for cimr jobs.')
    outfile = str(outdir) + '/' + str(args.out) + '_'
    logging.info(f' file prefix {str(args.out)} will be used for output files')

    data_type = args.data_type

    if args.process:
        if data_type in data_types:
            if args.file_name is not None:
                outfile = outfile + data_type + '.txt'
                infile = Infiler(data_type, args.file_name, args.genome_build)
                infile.read_file()
                genes = list(infile.list_genes())
                if data_type == 'eqtl':
                    queried = Querier(genes)
                    queried.form_query()
                    print()
                infile.write_file(outfile)
            else:
                logging.error(f' no file_name provided. nothing to process.')

    elif args.query:
        if data_type == 'gene':
            with open(args.file_name) as f:
                genes = f.read().splitlines()
            queried = Querier(genes)
            queried.form_query()
            if args.write_json is not None:
                annot_gene_file = str(outdir) + '/' + str(args.write_json)
                queried.write_json(annot_gene_file)
            if args.write_gene is not None:
                annot_gene_file = str(outdir) + '/' + str(args.write_gene)
                queried.write_gene(annot_gene_file)
        elif data_type == 'tad':        
            annotate_tad(args.file_name)
        else:
            logging.error(f' data-type is not recognized.')

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
        logging.error(f' data_type or file_name is not recognized. nothing to do.')
    

