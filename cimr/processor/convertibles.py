#!/usr/bin/env python3
"""Utility functions to convert between values and units
for downstream analyses.
e.g. log(OR) -> effect_size
"""

__author__ = "yoson park"


import re
import numpy
import pandas
import logging

from scipy import stats

from .constants import EFFECT_SIZE
from .constants import ODDS_RATIO
from .constants import ZSCORE
from .constants import PVALUE
from .constants import STANDARD_ERROR
from .constants import EFFECT_DIRECTION

from ..defaults import VERY_SMALL_P


def get_effect_direction(data):
    """Distinguish datasets with absolute beta effects only.
    data is assumed to be a pandas dataframe
    """
    effect_direction = None

    if EFFECT_SIZE in data:
        logging.debug(f' checking direction of effects.')
        effect_direction = numpy.sign(data[EFFECT_SIZE])
    elif EFFECT_DIRECTION in data:
        logging.debug(f' checking effect_dir column.')
        effect_direction = data[EFFECT_DIRECTION]
        # Make it equivalent to numpy.sign
        effect_direction = effect_direction.apply(
            lambda x: 1.0 if (x == '+' or x==1.0) else -1.0
        )

    if effect_direction is None:
        logging.warning(f' effect direction is not provided.')
    return effect_direction


def get_z(data, VERY_SMALL_P):
    """Given dataset with effect_size and standard_error or
    pvalue columns, calculate zscore.
    """
    if ZSCORE in data:
        logging.info(f' data contains zscore column.')
        return data

    else:
        z = None
        if PVALUE in data:
            logging.info(f' calculating zscore from pvalue.')
            z = convert_p_to_z(data, VERY_SMALL_P)
        elif STANDARD_ERROR in data and EFFECT_SIZE in data:
            logging.info('calculating zscore from se and beta')
            z = data[EFFECT_SIZE] / data[STANDARD_ERROR]

    if z is None:
        logging.warning(f' zscore could not be calculated from available data.')
        data[ZSCORE] = z
        return data


def convert_p_to_z(data, VERY_SMALL_P):
    """Calculate zscores from pvalues"""
    p = data[PVALUE].values

    if numpy.any(p == 0):
        logging.warning(f' pvalue column contains zero(s). This may be caused by numerical resolution limits. Consider using beta/se columns or check your input data.')

    effect_direction = get_effect_direction(data)
    abs_z = -stats.norm.ppf(p / 2)

    if numpy.any(numpy.isinf(abs_z)) and VERY_SMALL_P:
        logging.warning(' thresholding zscores.')
        min_p = numpy.min(p[numpy.logical_and(numpy.isfinite(abs_z), p != 0)])
        if VERY_SMALL_P < min_p:
            min_p = VERY_SMALL_P
        fix_z = -stats.norm.ppf(min_p / 2)
        logging.warning(f' using {fix_z} to fill in divergent zscores.')
        abs_z[numpy.isinf(abs_z)] = fix_z

    z = abs_z * effect_direction
    return z


def convert_z_to_p(zscore):
    """Given zscore, calculate pvalue."""
    return 2 * stats.norm.sf(numpy.absolute(zscore))


def convert_or_to_beta(odd_ratio):
    """Checks odds_ratio column values and returns beta."""
    if numpy.any(numpy.where(odd_ratio < 0)):
        logging.error(f' odds_ratio column includes negative values.')
    if numpy.any(numpy.where(odd_ratio == 0)):
        logging.warning(f' odds_ratio column includes zeroes.')
    return numpy.log(odd_ratio)


def estimate_se(effect_size, zscore):
    """Given effect_size and zscore, calculate standard_error."""
    return effect_size / zscore

