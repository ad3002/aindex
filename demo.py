#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.01.2018
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

from aindex import load_aindex

settings = {
  "index_prefix": "tests/kmers.23",
  "aindex_prefix": "tests/kmers.23",
  "reads_file": "tests/reads.reads",
}

index = load_aindex(settings)


