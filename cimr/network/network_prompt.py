#!/usr/bin/env python3

import logging
import network_similarity

def network_cli(args):
    args_dict = vars(args)
    outdir = args_dict.outdir 
    outdir.mkdir(exist_ok=True)
    if args_dict.random:
        if args_dict.randomselect is not None:
            randomselect = args_dict.randomselect
            network_similarity.write_subset(celltype, filesize, randomselect)
    pass


if __name__ == '__main__':
