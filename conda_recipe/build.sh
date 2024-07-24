#!/bin/bash

cd $SRC_DIR

make install

$PYTHON -m pip install --upgrade pip

$PYTHON -m pip install build installer

$PYTHON -m build

$PYTHON -m installer dist/aindex-{{ version }}-py3-none-any.whl
