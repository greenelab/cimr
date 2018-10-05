#!/usr/bin/env python

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
import logging
import warnings
import importlib

import cimr

class Twas():
    pass


class Netwas():
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
            default='output',
            help='path to directory where output files will be written to',
        )   
    subparsers.required = True
    subparsers.dest = 'subcommand'
    add_subparser_processor(subparsers)
    add_subparser_twas(subparsers)
    add_subparser_netwas(subparsers)
    args = parser.parse_args()
    return args

def add_subparser_processor(subparsers):
    parser = subparsers.add_parser(
        name='processor', help='process and integrate new data',
        description='process pull requests containing new association summary statistics '
                    'or new annotation files',
    )
    parser.add_argument(
        '--eqtl', 
        help='process association summary statistics file from expression- '
             'quantitative trait loci mapping',
    )
    parser.add_argument(
        '--gwas', 
        help='process association summary statistics file from genome-wide '
              'association studies mapping',
    )
    parser.add_argument(
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
    parser.add_argument(
        '--mr', action='store_true',
        help='transcriptome-wide association study using mendelian randomization',
    )
    parser.add_argument(
        '--abf', action='store_true',
        help='transcriptome-wide association study using approximate bayes factor',
    )
    parser.set_defaults(function='cimr.twas.twas_prompt.twas_cli')

def add_subparser_netwas(subparsers):
    parser = subparsers.add_parser(
        name='netwas', help='netwas using cimr data',
        description='run network association study tools '
                    'include options --svm and --rwr',        
    )
    parser.add_argument(
        '--svm', action='store_true',
        help='network-association study analysis using support vector machines',
    )
    parser.add_argument(
        '--rwr', action='store_true',
        help='network-assication study analysis using random walk with restarts',
    )
    parser.set_defaults(function='cimr.netwas.netwas_prompt.netwas_cli')

def main():
    """main prompt of cimr"""
    args = parse_arguments()
    module_name, function_name = args.function.rsplit('.', 1)
    module = importlib.import_module(module_name)
    function = getattr(module, function_name)
    function(args)

if __name__ == '__main__':
    main()


