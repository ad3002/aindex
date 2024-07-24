#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.01.2018
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import aindex

prefix_path = "tests/raw_reads.101bp.IS350bp25"
kmer2tf = aindex.get_aindex(prefix_path)

kmer = "A"*23
rkmer = "T"*23
kid = kmer2tf.get_kid_by_kmer(kmer)
print(kmer2tf.get_kmer_info_by_kid(kid))
print(kmer2tf[kmer], kid, kmer2tf.get_kmer_by_kid(kid), len(kmer2tf.pos(kmer)), kmer2tf.get_strand(kmer), kmer2tf.get_strand(rkmer))
kmer = kmer2tf.get_read(0, 23, 0)
pos = kmer2tf.pos(kmer)[0]
print(pos)

print(kmer2tf.get_kid_by_kmer(kmer), kmer2tf.get_kid_by_kmer(rkmer))

print(kmer2tf.get_hash_size())

print(kmer2tf.get_read(0, 123, 0))

print(kmer2tf.get_read(0, 123, 1))


k = 23
for p in kmer2tf.pos(kmer):
  print(kmer2tf.get_read(p, p+k))
  
test_kmer = "TAAGTTATTATTTAGTTAATACT"
right_kmer = "AGTTAATACTTTTAACAATATTA"

print(kmer2tf[kmer])

sequence = kmer2tf.get_read(0, 1023, 0)

print("Task 1. Get kmer frequency")
for i, (kmer, tf) in enumerate(kmer2tf.iter_sequence_kmers(sequence)):
    print(f"Position {i} kmer {kmer} freq = {tf}")
  
print("Task 2. Iter read by read, print the first 20 reads")
for rid, read in kmer2tf.iter_reads():
    if rid == 20:
        break
    print(rid, read)

print("Task 3. Iter reads by kmer, returs (start, next_read_start, read, pos_if_uniq|None, all_poses)")
# for read in aindex.iter_reads_by_kmer(test_kmer, kmer2tf):
#     print(read)

print("Task 4. Get distances in reads for two kmers, returns a list of (rid, left_kmer_pos, right_kmer_pos) tuples.")
# print aindex.get_left_right_distances(test_kmer, right_kmer, kmer2tf)

print("Task 5. Get layout for kmer, returns (max_pos, reads, lefts, rights, rids, starts), for details see source code")
# max_pos, reads, lefts, rights, rids, starts = get_layout_for_kmer(right_kmer, kmer2tf)
print("Central layout:")
# for read in reads:
#     print(read)
print("Left flanks:")
# print(lefts)
print("Right flanks:")
# print(rights)

print("Task 6. Iter reads by sequence, returns (start, next_read_start, read, pos_if_uniq|None, all_poses)")
sequence = "AATATTATTAAGGTATTTAAAAAATACTATTATAGTATTTAACATA"
# for read in iter_reads_by_sequence(sequence, kmer2tf):
#     print(read)



