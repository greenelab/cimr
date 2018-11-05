#!/usr/bin/env python3

import logging
import pathlib

def process_cli(args):
    args_dict = vars(args)
    outdir = args_dict.outdir 
    outdir.mkdir(exist_ok=True)
    pass

