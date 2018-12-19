#!/usr/bin/env python3

"""utilities and common file checks used across different functions in
processor subunits"""

import sys
import pandas
import pathlib
import logging

def findfile(filename):
    """find the file indicated in the prompt"""
    if(pathlib.Path(filename).resolve()):
        filename = str(filename)
        return filename
    else:
        logging.error(f'no file {filename} found for processing.')
        exit()


