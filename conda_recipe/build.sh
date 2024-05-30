#!/bin/bash

# Ensure we start in the right directory
cd $SRC_DIR

# Run the Makefile
make install

# Install the Python package
$PYTHON setup.py install
