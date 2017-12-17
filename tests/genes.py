

import sys
sys.path.append("/mnt/guatemala/akomissarov/Boechera_spatifolia/aindex")

from aindex import *

settings = {
  "index_prefix": "/mnt/guatemala/akomissarov/Boechera_spatifolia/raw.23.L3",
  "aindex_prefix": "/mnt/guatemala/akomissarov/Boechera_spatifolia/raw.23.L3",
  "reads_file": "/mnt/guatemala/akomissarov/Boechera_spatifolia/raw.reads",
}

index = load_aindex(settings)