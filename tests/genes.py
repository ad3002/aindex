

import sys
sys.path.append("/mnt/guatemala/akomissarov/Boechera_spatifolia/aindex")

from aindex import *
from trseeker.seqio.fasta_file import sc_iter_fasta


if __name__ == '__main__':
    

    settings = {
      "index_prefix": "/mnt/guatemala/akomissarov/Boechera_spatifolia/raw.23.L3",
      "aindex_prefix": "/mnt/guatemala/akomissarov/Boechera_spatifolia/raw.23.L3",
      "reads_file": "/mnt/guatemala/akomissarov/Boechera_spatifolia/raw.reads",
      "gene_fasta": "/mnt/guatemala/akomissarov/Boechera_spatifolia/apr1.fa",
    }

    k = 23
    index = load_aindex(settings)

    for seq_obj in sc_iter_fasta(settings["gene_fasta"]):

        for i in xrange(seq_obj.length-k+1):
            kmer = seq_obj.sequence[i:i+k]

            print i, kmer, index[kmer]

            for read in iter_reads_by_kmer(kmer, index):
                print read

            raw_input("?")





