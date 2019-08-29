#!/usr/bin/env python3

"""Utilities and common file checks used across different
processor classes

(c) YoSon Park
"""


import sys
import pandas
import pathlib
import logging
import subprocess

from pandas.api.types import is_numeric_dtype

from ..defaults import MAXCHROM


def set_chrom_dict():
    """Make a dictionary to standardize chromosome IDs in input files."""
    chrom_dict = {
        str(i):'chr' + str(i) for i in range(1, MAXCHROM)
    }
    chrom_dict.update({
        'X':'chr23',
        'Y':'chr24',
        'M':'chr25',
        'MT':'chr25',
        'chrX':'chr23',
        'chrY':'chr24',
        'chrM':'chr25',
        'chrMT':'chr25'
    })
    return chrom_dict, MAXCHROM


def find_file(file_name):
    """Check if a file exists and exit if not."""
    if (pathlib.Path(file_name).resolve()):
        file_name = str(file_name)
        logging.info(f' found {file_name}.')
        return file_name
    else:
        logging.error(f' no file {file_name} found for processing.')
        sys.exit()


def check_numeric(data, col):
    """Check for numeric columns"""
    from pandas.api.types import is_numeric_dtype
    try:
        if is_numeric_dtype(data[col]):
            logging.info(f' {col} is numeric.')
        else:
            numdata = (data
                        .drop([col], axis=1)
                        .join(data[col].apply(pandas.to_numeric, errors='coerce'))
                        )
            numcol = numdata[col].isnull().values().sum()
            logging.warning(f' %s rows in %s are non-numeric' % (numcol, col,))
            logging.warning(f' {col} is tested by coercing into numeric values.')
            return numdata
    except:
        logging.error(f' the format of %s is not testable.' % (col,))
        print(data.head(n=2))
        sys.exit(1)


def intersect_set(list1, list2):
    """Make a list of intersect set values"""
    return (set(list1) & set(list2))


