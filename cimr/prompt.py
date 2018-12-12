#!/usr/bin/env python3

"""provides general options for summary-based medelian randomization 
and other transcriptome-wide association study (twas) methods.
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

class Twas():
    pass


class Network():
    pass


def parse_arguments():
    """parse command line arguments for subprocesses of cimr."""
    parser = argparse.ArgumentParser(
        description='cimr: continuous integration for twas')
    parser.add_argument('--version', action='version', version=f'v{cimr.__version__}')
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='cimr subcommands:',
    )
    for subparser in subparsers.choices.values():
        subparser.add_argument(
            '--outdir',
            default='outdir',
            help='path to directory where output files will be written to',
        )
        subparser.add_argument(
            '--log',
            default='WARNING',
            choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
            help='logging level for stderr logging',
        )      
    subparsers.required = True
    subparsers.dest = 'subcommand'
    add_subparser_processor(subparsers)
    add_subparser_twas(subparsers)
    add_subparser_network(subparsers)
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
        '--eqtl', 
        help='process association summary statistics file from expression- '
             'quantitative trait loci mapping',
    )
    pargs.add_argument(
        '--gwas', 
        help='process association summary statistics file from genome-wide '
              'association studies mapping',
    )
    pargs.add_argument(
        '--tad', 
        help='process annotations for topologically associated domains',
    )
    parser.set_defaults(function='cimr.processor.processor_prompt.process_cli')

def add_subparser_twas(subparsers):
    parser = subparsers.add_parser(
        name='twas', help='twas using cimr data',
        description='run gene-based analyses or transcriptome-wide association analysis '
                    'before network-wide association studies'
    )
    targs = parser.add_mutually_exclusive_group()
    targs.add_argument(
        '--mr', default=False,
        action='store_true',
        help='transcriptome-wide association study using 2-sample-based mendelian randomization',
    )
    targs.add_argument(
        '--abf', default=False,
        action='store_true',
        help='transcriptome-wide association study using approximate bayes factor',
    )
    parser.set_defaults(function='cimr.twas.twas_prompt.twas_cli')

def add_subparser_network(subparsers):
    parser = subparsers.add_parser(
        name='network', help='network analysis using cimr data',
        description='run network analysis tools '
                    'include options --svm and --rwr',        
    )
    nargs = parser.add_mutually_exclusive_group()
    nargs.add_argument(
        '--random', default=False,
        action='store_true',
        help='select random edges from a network. '
             'use --randomcount to indicate number of edges to select',
    )
    parser.add_argument(
        '--randomcount', default=100000,
        dest='randomcount',
        nargs=1,
        type=int,
        help='select indicated number of random edges from a network'
             'use when --random is selected',
    )
    parser.add_argument(
        '--celltype',
        dest='celltype',
        default='global',
        nargs=1,
        type=str,
        help='cell or tissue type that both represents the specificity '
             'of the given network and the file name'
             'e.g. for giant2, file name is assumed to be celltype.dat',
    )
    parser.add_argument(
        '--filesize',
        dest='filesize',
        default=10000000,
        nargs=1,
        type=int,
        help='size of the file containing the network in the format of '
             'edge0 edge1 weight',
    )
    nargs.add_argument(
        '--svm', default=False,
        action='store_true',
        help='network analysis using support vector machines',
    )
    nargs.add_argument(
        '--rwr', default=False,
        action='store_true',
        help='network analysis using random walk with restarts',
    )
    parser.set_defaults(function='cimr.network.network_prompt.network_cli')

def main():
    """main prompt of cimr"""
    args = parse_arguments()
    module_name, function_name = args.function.rsplit('.', 1)
    module = importlib.import_module(module_name)
    function = getattr(module, function_name)
    function(args)

if __name__ == '__main__':
    main()


