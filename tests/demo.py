#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.01.2018
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import os
from aindex import *

settings = {
  "index_prefix": "tests/test.23",
  "aindex_prefix": "tests/test.23",
  "reads_file": "tests/test.reads",
  "sdat_file": "tests/test.23.sdat",
}

script_folder = os.path.dirname(os.path.abspath(__file__))
command = f"python {script_folder}/compute_aindex.py -i {script_folder}/../tests/raw_reads.101bp.IS350bp25_1.fastq,{script_folder}/../tests/raw_reads.101bp.IS350bp25_2.fastq -t fastq -o {script_folder}/../tests/test --sort 1"
print(command)
os.system(command)

index = load_aindex(settings)

k = 23
sequence = "TAAGTTATTATTTAGTTAATACTTTTAACAATATTATTAAGGTATTTAAAAAATACTATTATAGTATTTAACATAGTTAAATACCTTCCTTAATACTGTTAAATTATATTCAATCAATACATATATAATATTATTAAAATACTTGATAAGTATTATTTAGATATTAGACAAATACTAATTTTATATTGCTTTAATACTTAATAAATACTACTTATGTATTAAGTAAATATTACTGTAATACTAATAACAATATTATTACAATATGCTAGAATAATATTGCTAGTATCAATAATTACTAATATAGTATTAGGAAAATACCATAATAATATTTCTACATAATACTAAGTTAATACTATGTGTAGAATAATAAATAATCAGATTAAAAAAATTTTATTTATCTGAAACATATTTAATCAATTGAACTGATTATTTTCAGCAGTAATAATTACATATGTACATAGTACATATGTAAAATATCATTAATTTCTGTTATATATAATAGTATCTATTTTAGAGAGTATTAATTATTACTATAATTAAGCATTTATGCTTAATTATAAGCTTTTTATGAACAAAATTATAGACATTTTAGTTCTTATAATAAATAATAGATATTAAAGAAAATAAAAAAATAGAAATAAATATCATAACCCTTGATAACCCAGAAATTAATACTTAATCAAAAATGAAAATATTAATTAATAAAAGTGAATTGAATAAAATTTTGAAAAAAATGAATAACGTTATTATTTCCAATAACAAAATAAAACCACATCATTCATATTTTTTAATAGAGGCAAAAGAAAAAGAAATAAACTTTTATGCTAACAATGAATACTTTTCTGTCAAATGTAATTTAAATAAAAATATTGATATTCTTGAACAAGGCTCCTTAATTGTTAAAGGAAAAATTTTTAACGATCTTATTAATGGCATAAAAGAAGAGATTATTACTATTCAAGAAAAAGATCAAACACTTTTGGTTAAAACAAAAAAAACAAGTATTAATTTAAACACAATTAATGTGAATGAATTTCCAAGAATAAGGTTTAATGAAAAAAACGATTTAAGTGAATTTAATCAATTCAAAATAAATTATTCACTTTTAGTAAAAGGCATTAAAAAAATTTTTCACTCAGTTTCAAATAATCGTGAAATATCTTCTAAATTTAATGGAGTAAATTTCAATGGATCCAATGGAAAAGAAATATTTTTAGAAGCTTCTGACACTTATAAACTATCTGTTTTTGAGATAAAGCAAGAAACAGAACCATTTGATTTCATTTTGGAGAGTAATTTACTTAGTTTCATTAATTCTTTTAATCCTGAAGAAGATAAATCTATTGTTTTTTATTACAGAAAAGATAATAAAGATAGCTTTAGTACAGAAATGTTGATTTCAATGGATAACTTTATGATTAGTTACACATCGGTTAATGAAAAATTTCCAGAGGTAAACTACTTTTTTGAATTTGAACCTGAAACTAAAATAGTTGTTCAAAAAAATGAATTAAAAGATGCACTTCAAAGAATTCAAACTTTGGCTCAAAATGAAAGAACTTTTTTATGCGATATGCAAATTAACAGTTCTGAATTAAAAATAAGAGCTATTGTTAATAATATCGGAAATTCTCTTGAGGAAATTTCTTGTCTTAAATTTGAAGGTTATAAACTTAATATTTCTTTTAACCCAAGTTCTCTATTAGATCACATAGAGTCTTTTGAATCAAATGAAATAAATTTTGATTTCCAAGGAAATAGTAAGTATTTTTTGATAACCTCTAAAAGTGAACCTGAACTTAAGCAAATATTGGTTCCTTCAAGATAATGAATCTTTACGATCTTTTAGAACTACCAACTACAGCATCAATAAAAGAAATAAAAATTGCTTATAAAAGATTAGCAAAGCGTTATCACCCTGATGTAAATAAATTAGGTTCGCAAACTTTTGTTGAAATTAATAATGCTTATTCAATATTAAGTGATCCTAACCAAAAGGAAAAATATGATTCAATGCTGAAAGTTAATGATTTTCAAAATCGCATCAAAAATTTAGATATTAGTGTTAGATGACATGAAAATTTCATGGAAGAACTCGAACTTCGTAAGAACTGAGAATTTGATTTTTTTTCATCTGATGAAGATTTCTTTTATTCTCCATTTACAAAAA"
test_kmer = "TAAGTTATTATTTAGTTAATACT"
right_kmer = "AGTTAATACTTTTAACAATATTA"

print("Task 1. Get kmer frequency")
# raw_input("\nReady?")
for i in range(len(sequence)-k+1):
    kmer = sequence[i:i+k]
    print("Position %s kmer %s freq = %s" % (i, kmer, index[kmer]))

print("Task 2. Iter read by read, print the first 20 reads")
# raw_input("\nReady?")
for i, read in enumerate(index.iter_reads()):
    if i == 20:
        break
    print(i, read)

print("Task 3. Iter overall reads by kmer, check the consistency")
with open(settings["sdat_file"]) as fh:
    for i, line in enumerate(fh):
        kmer, tf = line.strip().split("\t")
        tf = int(tf)
        print(kmer, tf, index[kmer])
        if i == 20:
            break

with open(settings["sdat_file"]) as fh:
    for line in fh:
        kmer, tf = line.strip().split("\t")
        tf = int(tf)
        # print(kmer, tf, index[kmer])
        assert tf == index[kmer]

print("Task 4. Iter reads by kmer, returs (start, next_read_start, read, pos_if_uniq|None, all_poses)")
# raw_input("\nReady?")
for read in iter_reads_by_kmer(test_kmer, index):
    print(read)


print("Task 5. Get distances in reads for two kmers, returns a list of (rid, left_kmer_pos, right_kmer_pos) tuples.")
# raw_input("\nReady?")
print(get_left_right_distances(test_kmer, right_kmer, index))


print("Task 6. Get layout for kmer, returns (max_pos, reads, lefts, rights, rids, starts), for details see source code")
# raw_input("\nReady?")
max_pos, reads, lefts, rights, rids, starts = get_layout_for_kmer(right_kmer, index)
print("Central layout:")
for read in reads:
    print(read)
print("Left flanks:")
print(lefts)
print("Right flanks:")
print(rights)

print("Task 7. Iter reads by sequence, returтs (start, next_read_start, read, pos_if_uniq|None, all_poses)")
# raw_input("\nReady?")
sequence = "AATATTATTAAGGTATTTAAAAAATACTATTATAGTATTTAACATA"
for read in iter_reads_by_sequence(sequence, index):
    print(read)


print("Task 8. Iter reads by kmer with reads as SE, returns (start, next_read_start, subread, kmere_pos, -1|0|1 for spring_pos, was_reversed, poses_in_read)")
# raw_input("\nReady?")
user_reads = set()
sequence = "AATATTATTAAGGTATTTAAAAAATACTATTATAGTATTTAACATA"
for rid, nextrid, read, pos, spring_pos, was_reversed, poses in get_reads_se_by_kmer(kmer, index, user_reads, k=23):
    print(rid, read, pos)


