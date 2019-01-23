#!/bin/bash
set -e
set -x

if [ $@ == "cimr" ]; then
  python3 setup.py build
  python3 setup.py install
fi

