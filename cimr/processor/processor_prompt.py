#!/usr/bin/env python3

"""General prompt for the cimr processor subprocess.

(c) YoSon Park
"""


import sys
import pathlib
import logging


from .utils import Infiler
from .utils import Integrator

from .query import Querier
from .query import Snpper

from .tad import Tadpole

from .yamler import Yamler
from .yamler import load_yaml

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

    try:
        yaml_file = pathlib.Path(yaml_file)
        yaml_file_path = yaml_file.resolve(strict=True)
        logging.info(f' processing {yaml_file_path}')
        yaml_data = load_yaml(yaml_file)
        y = Yamler(yaml_data)
        y.check_data_file()
        data_type = y.data_type
        file_name = pathlib.Path(y.outfile_path)
        outdir = pathlib.Path('processed_data/' + str(data_type))
        pathlib.Path('processed_data/').mkdir(exist_ok=True)
        out_path = y.outfile_path.replace('submitted', 'processed')
        out_path = out_path

        return data_type, file_name, outdir, out_path

    except FileNotFoundError:
        logging.info(f' {yaml_file_path} is not accesible')
        sys.exit(1)


def processor_cli(args):
    """cimr processor subprocess cli"""

    if args.process:

        if args.yaml_file:
            data_type, file_name, outdir, out_path = convert_yaml(args.yaml_file)
        else:
            data_type = args.data_type
            file_name = args.file_name
            outdir = args.outdir
            out_path = str(outdir) + '/' + str(args.out)
            
        if check_type(data_type):

            outdir.mkdir(exist_ok=True)
            logging.info(f' output directory: {str(outdir)}')

            if file_name:
                outfile = out_path + '.txt.gz'
                infile = Infiler(
                    data_type, 
                    file_name, 
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
                sys.exit(1)
            
            # elif data_type == 'snp':
            #     Snpper()

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
                grow_tadpoles(args)
            
            logging.info(f' finished processing {file_name}')

        else:
            sys.exit(1)

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
        sys.exit(1)
    

