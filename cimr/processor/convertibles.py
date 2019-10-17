#!/usr/bin/env python3
"""Utility functions to convert between values and units
for downstream analyses.
e.g. log(OR) -> effect_size
"""

__author__ = "yoson park"


import re
import sys
import numpy
import pandas
import logging

from scipy import stats

from .constants import (EFFECT_SIZE, ODDS_RATIO, ZSCORE,
    PVALUE, STANDARD_ERROR, EFFECT_DIRECTION)

from ..defaults import VERY_SMALL_P


def get_effect_direction(data):
    """Distinguish datasets with absolute beta effects only.
    data is assumed to be a pandas dataframe.
    """
    effect_direction = None

    if EFFECT_SIZE in data:
        logging.debug(f' checking direction of effects.')
        effect_direction = numpy.sign(data[EFFECT_SIZE])
        logging.debug(f' effect_direction[0:3]: {effect_direction[0:3]}')
    elif EFFECT_DIRECTION in data:
        logging.debug(f' checking effect_dir column.')
        effect_direction = data[EFFECT_DIRECTION]
        # Make it equivalent to numpy.sign
        effect_direction = effect_direction.apply(
            lambda x: 1.0 if (x == '+' or x == 1.0) else -1.0
        )
        logging.debug(f' effect_direction[0:3]: {effect_direction[0:3]}.')

    if effect_direction is None:
        logging.warning(f' effect direction is not provided.')

    return effect_direction


def _p_to_z(data, VERY_SMALL_P):
    """Calculate z-scores from p-values."""
    p = data[PVALUE].values

    if numpy.any(p == 0):
        logging.warning(f' p-value column contains zero(s).')
        logging.warning(f' This may be caused by numerical resolution limits.')
        logging.warning(f' Consider using beta/se columns or check your input data.')

    effect_direction = get_effect_direction(data)
    abs_z = -stats.norm.ppf(p / 2)

    if numpy.any(numpy.isinf(abs_z)) and VERY_SMALL_P:
        logging.warning(' thresholding zscores.')
        min_p = numpy.min(p[numpy.logical_and(numpy.isfinite(abs_z), p != 0)])
        if VERY_SMALL_P < min_p:
            min_p = VERY_SMALL_P
        fix_z = -stats.norm.ppf(min_p / 2)
        logging.warning(f' using {fix_z} to fill in divergent z-scores.')
        abs_z[numpy.isinf(abs_z)] = fix_z

    z = abs_z * effect_direction
    return z


def _or_to_beta(odd_ratio):
    """Checks odds_ratio column values and returns beta."""
    if numpy.any(numpy.where(odd_ratio < 0)):
        logging.error(f' odds_ratio column includes negative values.')
        sys.exit(1)
    if numpy.any(numpy.where(odd_ratio == 0)):
        logging.error(f' odds_ratio column includes zeroes.')
        sys.exit(1)

    return numpy.log(odd_ratio)


def estimate_se(data):
    """Given z-scores and effect sizes, estimate standard errors."""
    if STANDARD_ERROR in data.columns:
        logging.debug(f' {STANDARD_ERROR} is included.')
    elif (ZSCORE in data.columns and
        EFFECT_SIZE in data.columns and
        STANDARD_ERROR not in data.columns):
        logging.info(f' estimating {STANDARD_ERROR} from {EFFECT_SIZE} and {ZSCORE}.')
        data[STANDARD_ERROR] = data[EFFECT_SIZE] / data[ZSCORE]
        data['se_est_from_z_beta'] = 'True'
    else:
        logging.error(f' {STANDARD_ERROR} cannot be estimated based on available data.')
        sys.exit(1)

    return data


def get_z(data):
    """Given dataset with effect_size and standard_error or
    p-value columns, calculate z-scores.
    """
    if ZSCORE in data:
        logging.info(f' data contains z-score column.')
        z = data[ZSCORE]
    else:
        z = None
        if PVALUE in data:
            logging.info(f' calculating z-scores from p-values.')
            z = _p_to_z(data, VERY_SMALL_P)
        elif STANDARD_ERROR in data and EFFECT_SIZE in data:
            logging.info('calculating z-scores from se and beta')
            z = data[EFFECT_SIZE] / data[STANDARD_ERROR]

    if z is None:
        logging.warning(f' z-scores could not be calculated based on available data.')

    return data


def convert_z_to_p(data):
    """Given z-scores, calculate p-values."""
    if PVALUE in data.columns:
        logging.debug(f' {PVALUE} is included.')
    elif (ZSCORE in data.columns and
        PVALUE not in data.columns):
        data[PVALUE] = 2 * stats.norm.sf(numpy.absolute(data[ZSCORE]))
        data['p_est_from_z'] = 'True'
    else:
        logging.error(f' {PVALUE} cannot be estimated based on available data.')
        sys.exit(1)

    return data


def convert_p_to_z(data):
    if (PVALUE in data.columns and
        ZSCORE not in data.columns):
        data[ZSCORE] = _p_to_z(data, VERY_SMALL_P)
        data['z_est_from_p'] = 'True'
    else:
        logging.debug(f' {ZSCORE} is included.')

    return data


def convert_or_to_beta(data):
    """Given odds ratio, calculate effect size."""
    if EFFECT_SIZE in data.columns:
        logging.debug(f' {EFFECT_SIZE} is included.')
    elif (ODDS_RATIO in data.columns and
        EFFECT_SIZE not in data.columns):
        data[EFFECT_SIZE] = _or_to_beta(data[ODDS_RATIO])
    else:
        logging.warning(f' {EFFECT_SIZE} could not be estimated from available data.')

    return data

