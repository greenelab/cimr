#!/usr/bin/env python3

"""The main prompt of cimr.

Options and parameters are listed for cimr subprocesses.

(c) YoSon Park
"""

__author__ = "YoSon Park"
__copyright__ = "Copyright 2018, YoSon Park"
__credits__ = ["YoSon Park"]
__license__ = "BSD"
__maintainer__ = "YoSon Park"
__status__ = "production"

import os
import sys
import argparse
import pathlib
import logging
import warnings
import importlib

import cimr


def parse_arguments():
    """Parse command line arguments for subprocesses of cimr."""
    parser = argparse.ArgumentParser(
        description='cimr: continuous integration of summary statistics files for network analysis')
    parser.add_argument('--version', action='version', version=f'v{cimr.__version__}')
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
            '--outdir',
            default='outdir',
            dest='outdir',
            nargs='?',
            type=pathlib.Path,
            help='path to directory where output files will be written to',
        )
        subparser.add_argument(
            '--out',
            default='cimr',
            dest='out',
            nargs='?',
            type=pathlib.Path,
            help='prefix for output files',
        )
        subparser.add_argument(
            '--log',
            default='info',
            dest='loglevel',
            nargs='?',
            choices=['debug', 'info', 'warning', 'error', 'critical'],
            help='logging level for stderr logging.',
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
        '--process',
        default=False,
        dest='process',
        action='store_true',
        help='automated checks of the input data to be used for --integreate',
    )
    pargs.add_argument(
        '--integrate',
        default=False,
        dest='integrate',
        action='store_true',
        help='integration of the input data into the data repository',
    )
    pargs.add_argument(
        '--query',
        default=False,
        dest='query',
        action='store_true',
        help='pull down annotations for a list of gene_id\'s',
    )

    # common processor arguments
    parser.add_argument(
        '--file-name',
        default=None,
        type=pathlib.Path,
        dest='file_name',
        help='file containing summary statistics or annotation data',
    )
    parser.add_argument(
        '--data-type',
        default=None,
        dest='data_type',
        help='currently supported datatypes include \{\'gwas\', \'eqtl\'\}',
    )
    parser.add_argument(
        '--genome-build',
        default='b38',
        dest='genome_build',
        help='human genome build used for the input file mapping',
    )

    # query-specific arguments
    parser.add_argument(
        '--write-json',
        default=None,
        type=pathlib.Path,
        dest='write_json',
        help='write results of the gene annotation query as a json file.',
    )
    parser.add_argument(
        '--write-gene',
        default=None,
        type=pathlib.Path,
        dest='write_gene',
        help='write results of the gene annotation query as a text file. '
             'the output will include the following: '
             'official gene symbol, entrez gene id, and ensembl gene id.',
    )
    
    # integrate-specific arguments
    parser.add_argument(
        '--can-be-public',
        default=True,
        dest='can_be_public',
        help='a boolean variable indicating whether the integrated data can '
             'be made public. default is True.',
    )
    parser.add_argument(
        '--temp-dir',
        default='cimr-adb-temp',
        dest='temp_dir',
        help='temporary directory name to clone cimr database into.',
    )
    
    parser.set_defaults(function='cimr.processor.processor_prompt.processor_cli')


def add_subparser_gene(subparsers):
    parser = subparsers.add_parser(
        name='gene', help='wrapper for gene-based analysis',
        description='run gene-based analyses or transcriptome-wide association analysis '
                    'and generate gene-based scores'
    )
    targs = parser.add_mutually_exclusive_group()
    targs.add_argument(
        '--mr', 
        default=False,
        action='store_true',
        help='association study using two-sample-based mendelian randomization',
    )
    targs.add_argument(
        '--abf', 
        default=False,
        action='store_true',
        help='colocalization test using approximate bayes factor',
    )
    parser.set_defaults(function='cimr.gene.gene_prompt.gene_cli')


def add_subparser_network(subparsers):
    parser = subparsers.add_parser(
        name='network', help='network analysis using cimr data',
        description='run network analysis tools ',        
    )

    parser.add_argument(
        '--random-count', 
        default=100000,
        dest='random_count',
        nargs='?',
        type=int,
        help='select indicated number of random edges from a network'
             'use when --random is selected',
    )
    parser.add_argument(
        '--cell-type',
        dest='cell_type',
        default='global',
        nargs='?',
        type=str,
        help='cell or tissue type that both represents the specificity '
             'of the given network and the file name'
             'e.g. for giant2, file name is assumed to be celltype.dat',
    )
    parser.add_argument(
        '--filesize',
        dest='filesize',
        default=10000000,
        nargs='?',
        type=int,
        help='size of the file containing the network in the format of '
             'edge0 edge1 weight',
    )

    nargs = parser.add_mutually_exclusive_group()
    nargs.add_argument(
        '--random', 
        default=False,
        action='store_true',
        help='select random edges from a network. '
             'use --randomcount to indicate number of edges to select',
    )
    nargs.add_argument(
        '--svm', 
        default=False,
        action='store_true',
        help='network analysis using support vector machines',
    )
    nargs.add_argument(
        '--rwr', 
        default=False,
        action='store_true',
        help='network analysis using random walk with restarts',
    )

    parser.set_defaults(function='cimr.network.network_prompt.network_cli')


def main():
    """The main CLI prompt of cimr"""
    args = parse_arguments()
    loglevel = args.loglevel
    numeric_level = getattr(logging, loglevel.upper(), None)
    if loglevel in ['debug', 'info', 'warning', 'error', 'critical']:
        logging.basicConfig(level=numeric_level)
        if not isinstance(numeric_level, int):
            raise ValueError(' invalid log level: %s' % loglevel)
        logging.basicConfig(level=numeric_level)
    else:
        logging.error(f' --log argument must be debug, info, warning, error, or critical.')
        logging.error(f' --log level is set to \'info\' by default.')
    module_name, function_name = args.function.rsplit('.', 1)
    module = importlib.import_module(module_name)
    function = getattr(module, function_name)
    function(args)


if __name__ == '__main__':
    main()


