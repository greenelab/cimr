"""Update genomic positions to hg38 reference.
Uses pyliftover 0.4 currently, which only allows
lifting over point coordinates and not indels, etc.
"""


__author__ = 'yoson park'

import pandas
import logging
import pyliftover

from ..defaults import HG19TO38


def convert_coords(lifting, chrom, pos):
    """Using pyliftover to convert genomic positions from
    older reference to the GRCh38/hg38 reference.
    """
    _lifted_chrom = 'NA'
    _lifted_pos = 'NA'
    try:
        pos = int(pos)
        lifted = lifting.convert_coordinate(chrom, pos)
        if lifted:
            if len(lifted) > 1:
                logging.warning(f' Liftover with more than one candidate: {t.variant_id}')
            _lifted_chrom = lifted[0][0]
            _lifted_pos = lifted[0][1]
    except:
        pass
    return _lifted_chrom, _lifted_pos


def call_liftover(df):
    """Call pyliftover.LiftOver to update genomic coordinates."""
    logging.info(f' Updating genomic coordinates.')
    lifting = pyliftover.LiftOver(HG19TO38)
    new_chrom = []
    new_pos = []
    for t in df.itertuples():
        _lifted_chrom, _lifted_pos = convert_coords(lifting, t.chrom, t.pos)
        new_chrom.append(_lifted_chrom)
        new_pos.append(_lifted_pos)

    df = df.assign(chrom=new_chrom)
    df = df.assign(pos=new_pos)
    logging.info(f' str{df.shape[0]} variants after liftover')
    return df

