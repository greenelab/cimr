#!/usr/bin/env python3


import os
import sys
import pathlib
import logging

from cimr.processor.util import readfile


if __name__ == '__main__':
    filename = sys.argv[1]
    sumdata = readfile(filename)

