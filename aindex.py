#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.03.2015
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import os
import ctypes
from ctypes import cdll
from ctypes import *

from settings import dll_paths

for dll_path in dll_paths:
    if os.path.isfile(dll_path):
        print(dll_path, os.path.isfile(dll_path))
        lib = ctypes.CDLL(dll_path)
        break
else:
    raise Exception("Ariadna's dll was not found: %s" % str(dll_paths))

import numpy as np
import mmap
from collections import defaultdict

def get_revcomp(sequence):
    '''Return reverse complementary sequence.

    >>> complementary('AT CG')
    'CGAT'

    '''
    c = dict(zip('ATCGNatcgn~[]', 'TAGCNtagcn~]['))
    return ''.join(c.get(nucleotide, '') for nucleotide in reversed(sequence))


def hamming_distance(s1, s2):
    """ Get Hamming distance: the number of corresponding symbols that differs in given strings.
    """
    return sum(i != j for (i,j) in zip(s1, s2) if i != 'N' and j != 'N')


lib.AindexWrapper_get.argtypes = [c_void_p, c_char_p]
lib.AindexWrapper_get.restype = c_size_t

class AIndex(object):
    ''' Wrapper for working with cpp aindex implementation.
    '''

    local_buffer = {}

    def __init__(self, index_prefix):
        ''' Init Aindex wrapper and load perfect hash.
        '''
        self.obj = pointer(c_void_p(lib.AindexWrapper_new()))
        lib.AindexWrapper_load(self.obj, index_prefix)

    def __getitem__(self, kmer):
        ''' Return tf for kmer.
        '''
        return lib.AindexWrapper_get(self.obj, kmer)

    def get_kid_by_kmer(self, kmer):
        return lib.AindexWrapper_get_kid_by_kmer(self.obj, kmer)

    def get_kmer_by_kid(self, kid):
        lib.AindexWrapper_get_kmer_by_kid.restype = c_char_p
        return lib.AindexWrapper_get_kmer_by_kid(self.obj, kid)

    def load(self, index_prefix, max_tf):
        ''' Load aindex. max_tf limits 
        the size of returning array with positions.
        '''
        print("Loadind aindex: %s.*" % index_prefix)
        self.max_tf = max_tf
        lib.AindexWrapper_load_index(self.obj, c_char_p(index_prefix), c_int64(max_tf))

    def load_reads(self, reads_file):
        ''' Load reads with mmap and with aindex.
        '''
        print("Loadind reads with mmap: %s" % reads_file)
        with open(reads_file, "r+b") as f:
            self.reads = mmap.mmap(f.fileno(), 0)
            self.reads_size = self.reads.size()
        lib.AindexWrapper_load_reads(self.obj, reads_file)
        print("\tloaded %s chars." % self.reads_size)
    
    def iter_reads(self):
        ''' Iter over reads 
        and yield (start_pos, next_read_pos, read).
        '''
        start = 0
        end = 0
        N = len(self.reads)
        while True:
            while self.reads[end] != "\n":
                end += 1
            yield start, end+1, self.reads[start:end]
            start = end+1
            end = end+1
            if end >= N:
                break

    def iter_reads_se(self):
        ''' Iter over reads 
        and yield (start_pos, next_read_pos, 0|1|..., read).
        '''
        start = 0
        end = 0
        N = len(self.reads)
        rid = 0
        while True:
            while self.reads[end] != "\n":
                end += 1
            splited_reads = self.reads[start:end].split("~")
            for i, subread in enumerate(splited_reads):
                yield rid, start, i, subread
            rid += 1
            start = end+1
            end = end+1
            if end >= N:
                break

    def get_hash_size(self):
        ''' Get hash size.
        ''' 
        return lib.AindexWrapper_get_hash_size(self.obj)

    def get_rid(self, pos):
        ''' Get read id by positions in read file.
        '''
        lib.AindexWrapper_get_rid.restype = c_size_t
        return c_size_t(lib.AindexWrapper_get_rid(self.obj, c_size_t(pos))).value

    def get_kmer(self, pos, k=23):
        ''' Get kmer, revcomp kmer and corresondent tf 
        for given position in read file.
        '''
        kmer = ctypes.c_char_p("N"*k)
        rkmer = ctypes.c_char_p("N"*k)
        tf = lib.AindexWrapper_get_kmer(self.obj, pos, pointer(kmer), pointer(rkmer))
        return kmer.value, rkmer.value, tf

    def pos(self, kmer):
        ''' Return array of positions for given kmer.
        '''
        n = self.max_tf
        r = (ctypes.c_size_t*n)()
        lib.AindexWrapper_get_positions(self.obj, r, kmer)
        poses = []
        ids = []
        poses_array = []
        for i in xrange(n):
            if r[i] > 0:
                poses_array.append(r[i]-1)
            else:
                break
        return poses_array

    def set(self, poses_array, kmer, batch_size):
        ''' Update kmer batch in case of fixed batches.
        '''
        print("WARNING: called a function with the fixed batch size.")
        n = batch_size*2
        r = (ctypes.c_size_t*n)()
        for i, (rid,pos) in enumerate(poses_array):
            r[i+batch_size] = ctypes.c_size_t(rid)
            r[i] = ctypes.c_size_t(pos)
        lib.AindexWrapper_set_positions(self.obj, r, kmer)


def load_aindex(settings, prefix=None, reads=None, aindex_prefix=None, skip_reads=False, skip_aindex=False):
    ''' Wrapper over aindex loading.
    Load:
    1. basic aindex with tf only;
    2. reads (if not skip_reads set);
    3. aindex (if not skip_aindex set);
    '''
    if "aindex_prefix" in settings and settings["aindex_prefix"] is None:
        skip_aindex = True
    if "reads_file" in settings and settings["reads_file"] is None:
        skip_reads = True

    if prefix is None:
        prefix = settings["index_prefix"]
    if reads is None and not skip_reads:
        reads = settings["reads_file"]

    if not "max_tf" in settings:
        print("default max_tf is 10000")
        settings["max_tf"] = 10000

    if aindex_prefix is None and not skip_aindex:
        aindex_prefix = settings["aindex_prefix"]
    
    kmer2tf = AIndex(prefix)
    kmer2tf.max_tf = settings["max_tf"]
    if not skip_reads:
        kmer2tf.load_reads(reads)
    if not skip_aindex:
        settings["max_tf"] = kmer2tf.load(aindex_prefix, settings["max_tf"])
    return kmer2tf


def get_rid2poses(kmer, kmer2tf):
    ''' Wrapper that handle case when two kmer hits in one read.
    Return rid->poses_in_read dictionary for given kmer. 
    In this case rid is the start position in reads file.
    '''
    poses = kmer2tf.pos(kmer)
    hits = defaultdict(list)
    for pos in poses:
        start = kmer2tf.get_rid(pos)
        hits[start].append(c_size_t(pos).value - start)
    return hits


def iter_reads_by_kmer(kmer, kmer2tf, used_reads=None, only_left=False, skip_multiple=True, k=23):
    ''' Yield 
        (start, next_read_start, read, pos_if_uniq|None, all_poses)

    - only_left: return only left reads
    - skip_multiple: skip if more then one hit in read

    '''
    rid2poses = get_rid2poses(kmer, kmer2tf)
    for rid in rid2poses:
        if used_reads and rid in used_reads:
            continue
        poses = rid2poses[rid]
        if skip_multiple:
            if len(poses) > 1:
                continue

        end = rid
        while True:
            if kmer2tf.reads[end] == '\n':
                break
            end += 1
        read = kmer2tf.reads[rid:end]
        
        pos = poses[0]
        is_multiple_hit = len(poses) > 1
        if read[pos:pos+k] != kmer:
            read = get_revcomp(read)
            poses = [len(read) - x - k for x in poses]
            ori_pos = pos
            pos = poses[0]
            assert read[pos:pos+k] == kmer
                
        if only_left:
            spring_pos = read.find("~")
            poses = [x for x in poses if x < spring_pos]
            if len(poses) == 1:
                yield [rid, end+1, read, poses[0], poses]
            elif not poses:
                continue
            else:
                yield [rid, end+1, read, None, poses]
        else:
            one_hit = None
            if len(poses) == 1:
                one_hit = poses[0]
            yield [rid, end+1, read, one_hit, poses]


def iter_reads_by_sequence(sequence, kmer2tf, used_reads=None, only_left=False, skip_multiple=True, k=23):
    ''' Yield reads containing sequence
        (start, next_read_start, read, pos_if_uniq|None, all_poses)

    TODO: more effective implementation than if sequence in read
    '''
    if len(sequence) >= k:
        kmer = sequence[:k]
        for data in iter_reads_by_kmer(kmer, kmer2tf, used_reads=used_reads, only_left=only_left, skip_multiple=skip_multiple, k=k):
            all_poses = data[-1]
            read = data[2]
            for pos in all_poses:
                if sequence in read:
                    yield data            
    else:
        yield None


def iter_reads_by_sequence_and_hamming(sequence, hd, kmer2tf, used_reads=None, only_left=False, skip_multiple=True, k=23):
    ''' Yield reads containing sequence
        (start, next_read_start, read, pos_if_uniq|None, all_poses)

    TODO: more effective implementation than if sequence in read
    '''
    if len(sequence) >= k:
        kmer = sequence[:k]
        n = len(sequence)
        for data in iter_reads_by_kmer(kmer, kmer2tf, used_reads=used_reads, only_left=only_left, skip_multiple=skip_multiple, k=k):
            all_poses = data[-1]
            read = data[2]
            for pos in all_poses:
                if len(read[pos:]) == n:
                    if hamming_distance(read[pos:], sequence) <= hd:
                        yield data            
    else:
        yield None


def get_reads_se_by_kmer(kmer, kmer2tf, used_reads, k=23):
    ''' Split springs and return subreads.

    The used_reads is a set of (start,spring_pos_type) tuples.

    The spring_pos is equal to is_right in case of PE data.

    Return list of:
    (start, next_read_start, subread, kmere_pos, -1|0|1 for spring_pos, was_reversed, poses_in_read)

    '''
    result = []
    hits = get_rid2poses(kmer, kmer2tf)
    rkmer = get_revcomp(kmer)

    for hit in hits:
        end = hit
        while True:
            if kmer2tf.reads[end] == '\n':
                break
            end += 1
        poses = hits[hit]
        read = kmer2tf.reads[hit:end]
        was_reversed = 0

        pos = poses[0]
        if read[pos:pos+k] != kmer:
            read = get_revcomp(read)
            poses = [len(read) - x - k for x in poses]
            pos = poses[0]
            was_reversed = 1
            if read[pos:pos+k] != kmer:
                print("Critical error kmer and ref are not equal:")
                print(read[pos:pos+k])
                print(kmer)
                continue
                
        spring_pos = read.find("~")

        if spring_pos == -1:
            result.append([hit, end+1, read, pos, -1, was_reversed, poses])
            continue

        left_poses = [x for x in poses if x < spring_pos]
        right_poses = [x-spring_pos-1 for x in poses if x > spring_pos]

        if len(left_poses) == 1:
            pos = left_poses[0]
            if (hit,0) in used_reads:
                continue
            result.append([hit, end+1, read[:spring_pos], pos, 0, was_reversed, left_poses])
        if len(right_poses) == 1:
            pos = right_poses[0]
            if (hit,1) in used_reads:
                continue
            result.append([hit, end+1, read[spring_pos+1:], pos, 1, was_reversed, right_poses])
    return result



def get_left_right_distances(left_kmer, right_kmer, kmer2tf):
    '''
    Return a list of (rid, left_kmer_pos, right_kmer_pos) tuples.

    TODO: Handle cases: 1. distance in one subread; 2. distance in pair.
    TODO: implement it more efficently.
    '''
    hits = defaultdict(list)
    for pos in kmer2tf.pos(left_kmer):
        start = kmer2tf.get_rid(pos)
        hits[start].append((0,pos-start))

    for pos in kmer2tf.pos(right_kmer):
        start = kmer2tf.get_rid(pos)
        hits[start].append((1,pos-start))

    results = []
    for rid, hit in hits.items():
        if len(hit) == 1:
            continue
        both_kmers_found = set([x[0] for x in hit])
        if both_kmers_found == 1:
            continue
        if len(hit) == 2:
            if hit[0][0] == 0:
                results.append((rid, hit[0][1], hit[1][1]))
            else:
                results.append((rid, hit[1][1], hit[0][1]))
        for left_pos in [x[1] for x in hit if x[0] == 0]:
            for right_pos in [x[1] for x in hit if x[0] == 1]:
                results.append((rid, left_pos, right_pos))
    return results


def get_layout_for_kmer(kmer, kmer2tf, used_reads=None, k=23):
    ''' Get flanked layout and left and right reads, or empty string if no read.

    - skip rids from used_reads

    seen_rids - track multiple hits from one spring.

    Return:
        - kmer start in reads
        - center reads as layout
        - left reads
        - right reads
        - rids list
        - starts list
    Or inline:
        (start_pos, center_layout, lefts, rights, rids, starts)

    '''
    max_pos = 0
    reads = []
    if used_reads is None:
        used_reads= set()
    seen_rids = set()
    lefts = []
    rights = []
    rids = []
    starts = []
    for (rid,nrid,read,pos,poses) in iter_reads_by_kmer(kmer, kmer2tf, used_reads, only_left=False, skip_multiple=False, k=k):
        if rid in seen_rids:
            continue
        seen_rids.add(rid)
        pos = poses[0]
        spring_pos = read.find("~")
        left,right = read.split("~")
        if pos < spring_pos:
            lefts.append("")
            rights.append(right)
            read = left
        else:
            lefts.append("")
            rights.append(left)
            pos = pos - len(left) - 1
            read = right
        max_pos = max(max_pos,pos)
        reads.append(read)
        starts.append(pos)
        rids.append(rid)
    max_length = max([len(x)+max_pos-starts[i] for i,x in enumerate(reads)])
    for i,read in enumerate(reads):
        reads[i] = 'N'*(max_pos-starts[i]) + read + "N" * (max_length-max_pos+starts[i]-len(read))
    return max_pos, reads, lefts, rights, rids, start_pos

