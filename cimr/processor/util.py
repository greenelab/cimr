#!/usr/bin/env python3

"""utilities and common file checks used across different functions in
processor subunits"""

import sys
import pandas
import pathlib
import logging

def findfile(filename):
    if(pathlib.Path(filename).resolve()):
        filename = str(filename)
        return filename
    else:
        logging.error(f'no file {filename} found for processing.\n')

def readfile(filename):
    # sumdata = pandas.read_csv(filename, sep='\t', header=0)
    # if not sumdata.empty:
    #     return sumdata
    # else:
        # logging.error(f'no content in uploaded file {filename}.\n')
    sumdata = open(filename).readlines()
    return sumdata

if __name__ == '__main__':
    filename = sys.argv[1]
    filename = findfile(filename)
    sumdata = readfile(filename)
    print(filename)
    print(sumdata)
    # print(sumdata.info())


