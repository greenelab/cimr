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
from .yamler import convert_yaml

from ..defaults import DATA_TYPES
from ..defaults import FILE_EXTENSION


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


def main(args, 
         data_type, 
         genome_build, 
         file_name, 
         outdir, 
         out_path, 
         columnset):
    """A standard set of functions for cimr -process argument
    
    Within cimr, main() may be called with a cimr call using
    a set of parameters or a yaml file.

    main() call Infiler class from cimr processor utils and
    runs standard processing steps for the Infiler object.
    """
    if check_type(data_type):

        if file_name:
            if not str(out_path).endswith(FILE_EXTENSION):
                outfile = pathlib.Path(str(out_path) + '.tsv.gz')
            else:
                outfile = pathlib.Path(out_path)

            infile = Infiler(
                data_type, 
                file_name, 
                genome_build, 
                args.update_rsid, 
                outfile,
                args.chunksize,
                columnset=columnset
            )
            
            infile.read_file()
            
            if data_type == 'eqtl':
                genes = list(infile.list_genes())
                queried = Querier(genes)
                queried.form_query()

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
            logging.info(f' output has been saved as {outfile}')

        else:
            logging.error(f' no file_name provided; nothing to process.')
            sys.exit(1)
        

def processor_cli(args):
    """cimr processor subprocess cli"""

    if args.process:

        if args.yaml_file:
            yaml_file = [pathlib.Path(args.yaml_file),]
            genome_build, fileset, columnset = convert_yaml(yaml_file)
            for _file in fileset:
                data_type = _file.split('/')[-2]
                file_name = _file.split('/')[-1]
                pathlib.Path('processed_data').mkdir(exist_ok=True)
                outdir = 'processed_data/' + str(data_type) + '/'
                logging.info(f' making dir: {outdir}')
                pathlib.Path(outdir).mkdir(exist_ok=True)
                out_path = outdir + file_name
                main(
                    args, 
                    data_type, 
                    genome_build, 
                    _file, 
                    outdir, 
                    out_path, 
                    columnset
                )

        else:
            data_type = args.data_type
            genome_build = args.genome_build
            file_name = args.file_name
            outdir = args.outdir
            pathlib.Path(outdir).mkdir(exist_ok=True)
            out_path = str(outdir) + '/' + str(args.out)
            columnset = {}
            main(
                args, 
                data_type, 
                genome_build, 
                file_name, 
                outdir, 
                out_path, 
                columnset
            )
        
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
        logging.error(f' data_type or file_name is not recognized; nothing to do.')
        sys.exit(1)
    

