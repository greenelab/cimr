"""Update genomic positions to hg38 reference.
Uses pyliftover 0.4 currently, which only allows
lifting over point coordinates and not indels, etc.
"""


__author__ = 'yoson park'

import os
import sys
import pandas
import logging
import pyliftover

from ..defaults import (HG19TO38, HG18TO38)


def convert_coords(lifting, t):
    """Using pyliftover to convert genomic positions from
    older reference to the GRCh38/hg38 reference.
    """
    _lifted_chrom = 'NA'
    _lifted_pos = 'NA'
    chrom = t.chrom
    pos = t.pos
    variant_id = t.variant_id
    try:
        pos = int(pos)
        lifted = lifting.convert_coordinate(chrom, pos)
        if lifted:
            if len(lifted) > 1:
                logging.warning(f' Liftover with more than one candidate: {variant_id}')
            _lifted_chrom = lifted[0][0]
            _lifted_pos = lifted[0][1]
    except:
        logging.debug(f' {variant_id} is not an integer.')

    return _lifted_chrom, _lifted_pos


def call_liftover(df):
    """Call pyliftover.LiftOver to update genomic coordinates."""
    logging.info(f' updating genomic coordinates.')
    build = df['build'][0]

    if (build == 'hg37') | (build == 'hg19') | (build == 'b37'):
        chain = HG19TO38
    elif (build == 'hg18') | (build == 'b18'):
        chain = HG18TO38
    else:
        logging.error(f' genome build information is not available.')
        sys.exit(1)

    lifting = pyliftover.LiftOver(chain)
    new_chrom = []
    new_pos = []
    df['chrom_'+build] = df['chrom']
    df['pos_'+build] = df['pos']
    df['variant_id_'+build] = df['variant_id']
    for t in df.itertuples():
        _lifted_chrom, _lifted_pos = convert_coords(lifting, t)
        new_chrom.append(_lifted_chrom)
        new_pos.append(_lifted_pos)

    df = df.assign(chrom=new_chrom)
    df = df.assign(pos=new_pos)

    # update build information in the dataframe
    df['build'] = 'b38'
    logging.info(f' {str(df.shape[0])} variants after liftover')

    return df
