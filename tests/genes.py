

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

    used_reads = set()
    for seq_obj in sc_iter_fasta(settings["gene_fasta"]):

        for i in xrange(seq_obj.length-k+1):
            kmer = seq_obj.sequence[i:i+k]
            tf = index[kmer]
            if not tf:
                continue

            print i, kmer, tf

            hits = []
            for data in get_reads_se_by_kmer(kmer, index, used_reads):
                start, next_read_start, subread, pos, spring_pos, was_reversed, poses_in_read = data
                used_reads.append((srart, spring_pos))
                hits.append((pos, subread, poses_in_read, was_reversed))

            max_pos = max([x[0] for x in hits])

            for pos, subread, poses_in_read, was_reversed in hits:
                print "N"*(max_pos-pos) + subread

            raw_input("?")





