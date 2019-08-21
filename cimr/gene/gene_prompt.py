#!/usr/bin/env python3

"""Argument parsing and function prompter for cimr gene subprocess."""


import logging
from .write_job import Writer


def gene_cli(args):
    """Main arguement parser for the gene subprocess."""

    outdir = args.outdir
    outdir = args.outdir
    outdir.mkdir(exist_ok=True)
    logging.info(f' directory {str(outdir)} will be used for cimr jobs.')
    outfile = str(outdir) + '/' + str(args.out)
    logging.info(f' file prefix {str(args.out)} will be used for output files')

    if args.pull_source is not None:
        pass
    elif args.job_dest is not None:
        writer_writes = Writer(
            args.cell_type,
            args.trait_name,
            args.feature_name,
            args.job_dest,
            args.outfile
        )
        writer_writes.get_tad_file_name()
    else:
        logging.error(f' define either -query or -write-job for cimr gene to run.')


