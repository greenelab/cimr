#!/usr/bin/env python3

"""General prompt for the cimr processor subprocess.

(c) YoSon Park
"""


import sys
import logging


from .utils import Infiler
from .utils import Integrator

from .query import Querier
from .query import Snpper

from .tad import Tadpole

from ..defaults import DATA_TYPES


def check_type(types):
    """Check data types"""
    if types in DATA_TYPES:
        return True
    else:
        logging.error(f' data type is not recognized.')
        return False


def grow_tadpoles(args):
    """Standard processing of a new tad annotation file."""
    tads = Tadpole(
        file_name=args.file_name,
        study_id=args.study_id,
        pub_id=args.pub_id, 
        species=args.species, 
        cell_type=args.cell_type,
        data_type=args.data_type,
        protocol=args.protocol
    )
    tads.read_file()
    tads.write_file()


def convert_yaml(yaml_file):
    """Convert yaml parameters to cimr arguments"""
    import pathlib
    from .yamler import Yamler
    from .yamler import load_yaml

    try:
        yaml_file = pathlib.Path(yaml_file)
        yaml_file_path = yaml_file.resolve(strict=True)
        logging.info(f' processing {yaml_file_path}')
        yaml_data = load_yaml(yaml_file)
        y = Yamler(yaml_data)
        y.check_data_file()
    except FileNotFoundError:
        logging.info(f' {yaml_file_path} is not accesible')
        sys.exit(1)


def processor_cli(args):
    """cimr processor subprocess cli"""

    outdir = args.outdir 
    outdir.mkdir(exist_ok=True)
    logging.info(f' output directory is {str(outdir)}')
    outfile = str(outdir) + '/' + str(args.out)
    logging.info(f' output file prefix is {str(args.out)}')

    data_type = args.data_type

    if args.process:
        
        if args.yaml_file:
            convert_yaml(args.yaml_file)

        if check_type(data_type):
            
            if args.file_name:
                outfile = outfile + '_' + data_type + '.txt.gz'
                infile = Infiler(
                    data_type, 
                    args.file_name, 
                    args.genome_build, 
                    args.update_rsid, 
                    outfile,
                    args.chunksize
                )
                infile.read_file()
                
                if data_type == 'eqtl':
                    genes = list(infile.list_genes())
                    queried = Querier(genes)
                    queried.form_query()
            
            else:
                logging.error(f' no file_name provided. nothing to process.')
                sys.exit(0)
        
        # elif data_type == 'snp':
        #     Snpper()

        elif data_type == 'gene':
            
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
            grow_tadpoles(args)
        else:
            pass

    elif args.integrate:
        if check_type(data_type):
            integrating = Integrator(
                data_type, 
                args.file_name, 
                can_be_public=args.can_be_public, 
                genome_build=args.genome_build
            )
            integrating.make_local_db(args.temp_dir)

    else:
        logging.error(f' data_type or file_name is not recognized. nothing to do.')
        sys.exit(0)
    

