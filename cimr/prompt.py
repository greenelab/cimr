#!/usr/bin/env python3

"""The main prompt of cimr.

Options and parameters are listed for cimr subprocesses.

(c) YoSon Park
"""


__author__ = "YoSon Park"
__copyright__ = "Copyright 2018-2019, YoSon Park"
__credits__ = ["YoSon Park"]
__license__ = "BSD"
__maintainer__ = "YoSon Park"
__status__ = "production"


import os
import sys
import json
import argparse
import pathlib
import logging
import warnings
import importlib

import cimr

from .defaults import CHUNKSIZE
from .defaults import DATA_TYPES



def parse_arguments():
    """Parse command line arguments for subprocesses of cimr."""
    parser = argparse.ArgumentParser(
        description='cimr: continuously integrated meta-resource'
    )
    parser.add_argument(
        '-version',
        action='version',
        version=f'v{cimr.__version__}'
    )
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='cimr subcommands:',
    )

    subparsers.required = True
    subparsers.dest = 'subcommand'
    add_subparser_processor(subparsers)
    add_subparser_gene(subparsers)
    add_subparser_network(subparsers)

    for subparser in subparsers.choices.values():
        subparser.add_argument(
            '-outdir',
            default='outdir',
            dest='outdir',
            nargs='?',
            type=pathlib.Path,
            help='path to directory where output files will be written to',
        )
        subparser.add_argument(
            '-out',
            default='cimr',
            dest='out',
            nargs='?',
            type=pathlib.Path,
            help='prefix for output files',
        )
        subparser.add_argument(
            '-log',
            default='info',
            dest='loglevel',
            nargs='?',
            choices=['debug', 'info', 'warning', 'error', 'critical'],
            help='logger arguement for stderr logging level.',
        )
    args = parser.parse_args()

    return args


def add_subparser_processor(subparsers):
    parser = subparsers.add_parser(
        name='processor', help='process and integrate new data',
        description='process pull requests containing new association summary statistics '
                    'or new annotation files',
    )

    pargs = parser.add_mutually_exclusive_group()
    pargs.add_argument(
        '-process',
        default=False,
        dest='process',
        action='store_true',
        help='automated checks of the input data to be used for -integrate',
    )
    pargs.add_argument(
        '-integrate',
        default=False,
        dest='integrate',
        action='store_true',
        help='integration of the input data into the data repository',
    )
    pargs.add_argument(
        '-query',
        default=False,
        dest='query',
        action='store_true',
        help='pull down annotations for a list of gene_id\'s',
    )

    # common processor arguments
    parser.add_argument(
        '-file-name',
        type=pathlib.Path,
        dest='file_name',
        help='file containing summary statistics or annotation data',
    )
    parser.add_argument(
        '-catalog-name',
        default='catalog.txt',
        type=pathlib.Path,
        dest='catalog_name',
        nargs='?',
        help='catalog file name to store metadata',
    )
    parser.add_argument(
        '-data-type',
        dest='data_type',
        help='input file data type',
        choices=DATA_TYPES,
    )
    parser.add_argument(
        '-genome-build',
        default='b38',
        dest='genome_build',
        help='human genome build used for the input file mapping',
    )
    parser.add_argument(
        '-update-map',
        default=False,
        dest='update_map',
        action='store_true',
        help='whether to update b37 data to b38 based on provided references',
    )
    parser.add_argument(
        '-update-rsid',
        default=False,
        dest='update_rsid',
        action='store_true',
        help='whether to update snp rs IDs based on refseq',
    )
    parser.add_argument(
        '-chunksize',
        default=CHUNKSIZE,
        dest='chunksize',
        type=int,
        help='number of rows in the input file to process at a time',
    )
    parser.add_argument(
        '-yaml-file',
        dest='yaml_file',
        help='yaml file containing information about data',
    )
    parser.add_argument(
        '-column-set',
        default='{"key":"value"}',
        dest='columnset',
        type=json.loads,
        nargs='?',
        help='dictionary containing corresponding header names, '
             'if different from the cimr default',
    )

    # not required for data-type gwas
    parser.add_argument(
        '-cell-type',
        default=None,
        dest='cell_type',
        help='cell or tissue type for the contributed data. '
             'including this info is recommended to integrate eQTL or TAD data.',
    )

    # query-specific arguments
    parser.add_argument(
        '-write-json',
        default=None,
        type=pathlib.Path,
        dest='write_json',
        help='write results of the gene annotation query as a json file.',
    )
    parser.add_argument(
        '-write-gene',
        default=None,
        type=pathlib.Path,
        dest='write_gene',
        help='write results of the gene annotation query as a text file. '
             'the output will include the following columns: '
             'official gene symbol, entrez gene id, and ensembl gene id.',
    )

    # integrate-specific arguments
    parser.add_argument(
        '-can-be-public',
        default=True,
        dest='can_be_public',
        help='a boolean variable indicating whether the integrated data can '
             'be made public. default is True.',
    )
    parser.add_argument(
        '-temp-dir',
        default='cimr-adb-temp',
        dest='temp_dir',
        help='temporary directory name to clone cimr database into.',
    )
    parser.add_argument(
        '-study-id',
        dest='study_id',
        help='unique identification number for the study representing the '
             'submitted data. If left blank, cimr will generate a random ID. '
             'An example includes a unique identifier provided by the gene '
             'expression omnibus (GEO).'
    )
    parser.add_argument(
        '-pub-id',
        dest='pub_id',
        help='DOI for the article describing the submitted data. '
    )
    parser.add_argument(
        '-species',
        dest='species',
        help='species. Default is homo_sapiens (also represented as 9606 or human).'
    )
    parser.add_argument(
        '-protocol',
        dest='protocol',
        help='name of the protocol used. For example, Hi-C for tad coordinates.'
    )

    parser.set_defaults(function='cimr.processor.processor_prompt.processor_cli')


def add_subparser_gene(subparsers):
    parser = subparsers.add_parser(
        name='gene', help='wrapper for gene-based analysis',
        description='provide scripts for gene-based analyses and transcriptome-wide '
                    'association analysis or provide results from pre-run analyses.',
    )
    targs = parser.add_mutually_exclusive_group()
    targs.add_argument(
        '-pull',
        dest='pull_source',
        help='pull down available gene-based association results.',
        choices=['twas'],
    )
    targs.add_argument(
        '-write-job',
        dest='job_dest',
        help='write example scripts for gene-based association studies.',
        choices=['annotation', 'enloc'],
    )

    parser.add_argument(
        '-trait',
        default='GLGC:LDL_Cholesterol',
        dest='trait_name',
        help='name of the trait for -query or -write-job. '
             'check documentation for the latest list of available traits '
             '(https://cimr.readthedocs.io/processing/data.html)',
        # '...' included to indicate there are more options available
        choices=['GLGC:LDL_Cholesterol', 'CARDIoGRAM_C4D:Coronary_Artery_Disease', '...']
    )
    parser.add_argument(
        '-cell-type',
        default='GTEx:Liver',
        dest='cell_type',
        help='name of the cell type for -query or -write-job. '
             'check documentation for the latest list of available tissues '
             '(https://cimr.readthedocs.io/processing/data.html)',
        # '...' included to indicate there are more options available
        choices=['GTEx:Liver:cis', 'GTEx:Whole_Blood:cis', 'TCGA:LIHC:cis', '...']
    )
    parser.add_argument(
        '-feature',
        default='gene:ENSG',
        dest='feature_name',
        help='feature space to run combined-variant analysis in. Default is '
             '\'gene\' for most methods. '
             'check documentation for the latest list of available feature space. '
             '(https://cimr.readthedocs.io/processing/data.html)',
        # '...' included to indicate there are more options available
        choices=['gene:ENSG', 'tad:Liver', 'tad:HepG2', 'ldblock', '...'],
    )

    parser.set_defaults(function='cimr.gene.gene_prompt.gene_cli')


def add_subparser_network(subparsers):
    parser = subparsers.add_parser(
        name='network', help='network analysis using cimr data',
        description='run network analysis tools ',
    )

    parser.add_argument(
        '-random-count',
        default=100000,
        dest='random_count',
        nargs='?',
        type=int,
        help='select indicated number of random edges from a network'
             'use when -random is selected',
    )
    parser.add_argument(
        '-cell-type',
        dest='cell_type',
        default='global',
        nargs='?',
        type=str,
        help='cell or tissue type that both represents the specificity '
             'of the given network and the file name'
             'e.g. for giant, file name is assumed to be \{celltype\}.dat',
    )
    parser.add_argument(
        '-file-size',
        dest='file-size',
        default=10000000,
        nargs='?',
        type=int,
        help='size of the file containing the network in the format of '
             'node1 node2 weight',
    )

    nargs = parser.add_mutually_exclusive_group()
    nargs.add_argument(
        '-random',
        default=False,
        action='store_true',
        help='select random edges from a network. '
             'use -randomcount to indicate number of edges to select',
    )
    nargs.add_argument(
        '-svm',
        default=False,
        action='store_true',
        help='network analysis using support vector machines',
    )

    parser.set_defaults(function='cimr.network.network_prompt.network_cli')


def set_log(args):
    """Set loglevel and format"""
    loglevel = args.loglevel
    numeric_level = getattr(logging, loglevel.upper(), None)

    FORMAT = '[%(asctime)-15s %(name)s:%(levelname)-8s] %(message)s'
    DATEFMT = '%Y-%m-%d %H:%M:%S'

    if loglevel in ['debug', 'info', 'warning', 'error', 'critical']:

        if not isinstance(numeric_level, int):
            raise ValueError(' invalid log level: %s' % loglevel)

        logging.basicConfig(
            level=numeric_level,
            format=FORMAT,
            datefmt=DATEFMT
        )
    else:
        raise ValueError(' -log argument must be debug, info, warning, error, or critical.')


def main():
    """The main CLI prompt of cimr"""
    args = parse_arguments()
    set_log(args)
    module_name, function_name = args.function.rsplit('.', 1)
    module = importlib.import_module(module_name)
    function = getattr(module, function_name)
    function(args)


if __name__ == '__main__':
    main()


