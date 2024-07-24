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
from intervaltree import IntervalTree
import mmap
from collections import defaultdict
import importlib.resources as pkg_resources

with pkg_resources.path('aindex.core', 'python_wrapper.so') as dll_path:
    dll_path = str(dll_path)

if not os.path.exists(dll_path):
    raise Exception(f"aindex's dll was not found: {dll_path}")

lib = cdll.LoadLibrary(dll_path)

lib.AindexWrapper_new.argtypes = []
lib.AindexWrapper_new.restype = c_void_p

lib.AindexWrapper_load.argtypes = [c_void_p, c_char_p, c_char_p]
lib.AindexWrapper_load.restype = None

lib.AindexWrapper_get.argtypes = [c_void_p, c_char_p]
lib.AindexWrapper_get.restype = c_size_t

lib.AindexWrapper_get_kid_by_kmer.argtypes = [c_void_p, c_char_p]
lib.AindexWrapper_get_kid_by_kmer.restype = c_size_t

lib.AindexWrapper_get_kmer_by_kid.argtypes = [c_void_p, c_size_t, c_char_p]
lib.AindexWrapper_get_kmer_by_kid.restype = None

lib.AindexWrapper_load_index.argtypes = [c_void_p, c_char_p, c_uint32]
lib.AindexWrapper_load_index.restype = None

lib.AindexWrapper_load_reads.argtypes = [c_void_p, c_char_p]
lib.AindexWrapper_load_reads.restype = None

lib.AindexWrapper_load_reads_index.argtypes = [c_void_p, c_char_p]
lib.AindexWrapper_load_reads_index.restype = None

lib.AindexWrapper_get_hash_size.argtypes = [c_void_p]
lib.AindexWrapper_get_hash_size.restype = c_size_t

lib.AindexWrapper_get_reads_size.argtypes = [c_void_p]
lib.AindexWrapper_get_reads_size.restype = c_size_t

lib.AindexWrapper_get_read.argtypes = [c_size_t, c_size_t, c_uint]
lib.AindexWrapper_get_read.restype = c_char_p

lib.AindexWrapper_get_read_by_rid.argtypes = [c_size_t]
lib.AindexWrapper_get_read_by_rid.restype = c_char_p

lib.AindexWrapper_get_rid.argtypes = [c_size_t]
lib.AindexWrapper_get_rid.restype = c_size_t

lib.AindexWrapper_get_start.argtypes = [c_size_t]
lib.AindexWrapper_get_start.restype = c_size_t

lib.AindexWrapper_get_strand.argtypes = [c_void_p]
lib.AindexWrapper_get_strand.restype = c_size_t

lib.AindexWrapper_get_kmer.argtypes = [c_void_p, c_size_t, c_char_p, c_char_p]
lib.AindexWrapper_get_kmer.restype = c_size_t

lib.AindexWrapper_get_positions.argtypes = [c_void_p, c_void_p, c_char_p]
lib.AindexWrapper_get_positions.restype = None

lib.AindexWrapper_set_positions.argtypes = [c_void_p, c_void_p, c_char_p]
lib.AindexWrapper_set_positions.restype = None


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


class AIndex(object):
    ''' Wrapper for cpp aindex implementation.
    '''

    obj = None
    loaded_header = False
    loaded_intervals = False
    loaded_reads = False
    
    def __init__(self, index_prefix):
        ''' Init Aindex wrapper and load perfect hash.
        '''
        self.obj = lib.AindexWrapper_new()
        if not (os.path.isfile(index_prefix + ".pf") and os.path.isfile(index_prefix + ".tf.bin") and os.path.isfile(index_prefix + ".kmers.bin")):
            raise Exception("One of index files was not found: %s" % str(index_prefix))
        tf_file = index_prefix + ".tf.bin"
        lib.AindexWrapper_load(self.obj, index_prefix.encode('utf-8'), tf_file.encode('utf-8'))

    def load(self, index_prefix, max_tf):
        ''' Load aindex. max_tf limits 
        the size of returning array with positions.
        '''
        print("Loadind aindex: %s.*" % index_prefix)

        if not (os.path.isfile(index_prefix + ".pf") and os.path.isfile(index_prefix + ".tf.bin") and os.path.isfile(index_prefix + ".kmers.bin") and os.path.isfile(index_prefix + ".index.bin") and os.path.isfile(index_prefix + ".indices.bin") and os.path.isfile(index_prefix + ".pos.bin")):
            raise Exception("One of index files was not found: %s" % str(index_prefix))

        self.max_tf = max_tf

        tf_file = index_prefix + ".tf.bin"

        lib.AindexWrapper_load_index(self.obj, index_prefix.encode('utf-8'), c_uint32(max_tf), index_prefix.encode('utf-8'), tf_file.encode('utf-8'))

    def load_local_reads(self, reads_file):
        ''' Load reads with mmap and with aindex.
        '''
        print("Loading reads with mmap: %s" % reads_file)
        with open(reads_file, "r+b") as f:
            self.reads = mmap.mmap(f.fileno(), 0)
            self.reads_size = self.reads.size()
        self.loaded_reads = True

    def load_reads_index(self, index_file, header_file=None):
        print("Loading reads index: %s" % index_file)
        self.rid2start = {}
        self.IT = IntervalTree()
        self.chrm2start = {}
        self.headers = {}
        with open(index_file) as fh:
            for line in fh:
                rid, start, end = line.strip().split("\t")
                self.rid2start[int(rid)] = (int(start), int(end))
                self.IT.addi(int(start), int(end), int(rid))
        self.loaded_intervals = True

        if header_file:
            print("Loading headers: %s" % header_file)
            with open(header_file) as fh:
                for rid, line in enumerate(fh):
                    head, start, length = line.strip().split("\t")
                    start = int(start)
                    length = int(length)
                    self.headers[rid] = head
                    chrm = head.split()[0].split(".")[0]
                    self.chrm2start[chrm] = start
                    self.IT.addi(start, start+length, head)
            self.loaded_header = True

    def load_reads(self, reads_file):
        ''' Load reads with mmap and with aindex.
        '''
        if not os.path.isfile(reads_file):
            raise Exception("Reads files was not found: %s" % str(reads_file))

        print("Loading reads with mmap: %s" % reads_file)
        lib.AindexWrapper_load_reads(self.obj, reads_file.encode('utf-8'))
        self.reads_size = lib.AindexWrapper_get_reads_size(self.obj)
        print(f"\tloaded {self.reads_size} chars.")

    def get_hash_size(self):
        ''' Get hash size.
        ''' 
        return lib.AindexWrapper_get_hash_size(self.obj)


    ### Getters for kmers

    def __getitem__(self, kmer):
        ''' Return tf for kmer.
        '''
        return lib.AindexWrapper_get(self.obj, kmer.encode('utf-8'))

    def get_strand(self, kmer):
        ''' Return strand for kmer:
            1 - the same as given
            2 - reverse complement
            0 - not found
        '''
        return lib.AindexWrapper_get_strand(self.obj, kmer.encode('utf-8'))

    def get_kid_by_kmer(self, kmer):
        ''' Return kmer id for kmer
        '''
        return lib.AindexWrapper_get_kid_by_kmer(self.obj, kmer.encode('utf-8'))
    
    def get_kmer_by_kid(self, kid, k=23):
        ''' Return kmer by kmer id 
        '''
        s = "N"*k
        kmer = ctypes.c_char_p()
        kmer.value = s.encode("utf-8")
        lib.AindexWrapper_get_kmer_by_kid(self.obj, c_size_t(kid), kmer)
        return kmer.value

    def get_kmer_info_by_kid(self, kid, k=23):
        ''' Get kmer, revcomp kmer and corresondent tf 
        for given kmer id.
        '''

        s = "N"*k

        kmer = ctypes.c_char_p()
        kmer.value = s.encode("utf-8")

        rkmer = ctypes.c_char_p()
        rkmer.value = s.encode("utf-8")

        tf = lib.AindexWrapper_get_kmer(self.obj, kid, kmer, rkmer)
        return kmer.value, rkmer.value, tf

    ### Getters for reads

    def get_rid(self, pos):
        ''' Get read id by positions in read file.
        '''
        return c_size_t(lib.AindexWrapper_get_rid(self.obj, c_size_t(pos))).value
    
    def get_start(self, pos):
        ''' Get read id by positions in read file.
        '''
        return c_size_t(lib.AindexWrapper_get_start(self.obj, c_size_t(pos))).value
    
    def get_read_by_rid(self, rid):
        ''' Get read sequence as string by rid.
        '''
        return lib.AindexWrapper_get_read_by_rid(self.obj, rid).decode("utf-8")

    def get_read(self, start, end, revcomp=False):
        ''' Get read by start and end positions.
        '''
        return lib.AindexWrapper_get_read(self.obj, start, end, revcomp).decode("utf-8")
        
    def iter_reads(self):
        ''' Iter over reads 
        and yield (start_pos, next_read_pos, read).
        '''
        if self.reads_size == 0:
            raise Exception("Reads were not loaded.")
        
        for rid in range(self.reads_size):
            yield rid, self.get_read_by_rid(rid)

    def iter_reads_se(self):
        ''' Iter over reads 
        and yield (start_pos, next_read_pos, 0|1|..., read).
        '''
        if self.reads_size == 0:
            raise Exception("Reads were not loaded.")
        
        for rid in range(self.reads_size):
            read = self.get_read_by_rid(rid)
            splited_reads = read.split("~")
            for i, subread in enumerate(splited_reads):
                yield rid, i, subread          

    def pos(self, kmer):
        ''' Return array of positions for given kmer.
        '''
        n = self.max_tf
        r = (ctypes.c_size_t*n)()
        kmer = str(kmer)
        lib.AindexWrapper_get_positions(self.obj, pointer(r), kmer.encode('utf-8'))
        poses_array = []
        for i in range(n):
            if r[i] > 0:
                poses_array.append(r[i]-1)
            else:
                break
        return poses_array

    def get_header(self, pos):
        ''' Get header information for position.
        '''
        if not self.loaded_header:
            return None
        rid = list(self.IT[pos])[0][2]
        return self.headers[rid]
    
    def iter_sequence_kmers(self, sequence):
        ''' Iter over kmers in sequence.
        '''
        for i in range(len(sequence)-23+1):
            kmer = sequence[i:i+23]
            yield kmer, self[kmer]

    def get_sequence_coverage(self, seq, cutoff=0):
        '''
        '''
        coverage = [0] * len(seq)
        for i in range(len(seq)-self.k+1):
            kmer = seq[i:i+23]
            tf = self[kmer]
            if tf >= cutoff:
                coverage[i] = tf
        return coverage
                
    def print_sequence_coverage(self, seq, cutoff=0):
        '''
        '''
        for i, tf in enumerate(self.get_sequence_coverage(seq, cutoff)):
            kmer = seq[i:i+23]
            print(i, kmer, tf)

    def get_rid2poses(self, kmer):
        ''' Wrapper that handle case when two kmer hits in one read.
        Return rid->poses_in_read dictionary for given kmer. 
        In this case rid is the start position in reads file.
        '''
        poses = self.pos(kmer)
        hits = defaultdict(list)
        for pos in poses:
            start = self.get_rid(pos)
            hits[start].append(c_size_t(pos).value - start)
        return hits

    ### Aindex manipulation

    def set(self, poses_array, kmer, batch_size):
        ''' Update kmer batch in case of fixed batches.
        '''
        print("WARNING: called a function with the fixed batch size.")
        n = batch_size*2
        r = (ctypes.c_size_t*n)()
        for i, (rid,pos) in enumerate(poses_array):
            r[i+batch_size] = ctypes.c_size_t(rid)
            r[i] = ctypes.c_size_t(pos)

        lib.AindexWrapper_set_positions(self.obj, pointer(r), kmer.encode('utf-8'))


def get_aindex(prefix_path, skip_aindex=False, max_tf=1_000_000):
    settings = {
        "index_prefix": f"{prefix_path}.23",
        "aindex_prefix": f"{prefix_path}.23",
        "reads_file": f"{prefix_path}.reads",
        "max_tf": max_tf,
        }
    if not os.path.isfile(prefix_path + ".23.pf"):
      print("No file", prefix_path + ".23.pf")
      return None
    if not os.path.isfile(prefix_path + ".23.tf.bin"):
      print("No file", prefix_path + ".23.tf.bin")
      return None
    if not os.path.isfile(prefix_path + ".23.kmers.bin"):
      print("No file", prefix_path + ".23.kmers.bin")
      return None
    if not skip_aindex:
        if not os.path.isfile(prefix_path + ".23.index.bin"):
            print("No file", prefix_path + ".23.index.bin")
            return None
        if not os.path.isfile(prefix_path + ".23.indices.bin"):
            print("No file", prefix_path + ".23.indeces.bin")
            return None
        if not os.path.isfile(prefix_path + ".23.pos.bin"):
            print("No file", prefix_path + ".23.pos.bin")
            return None
        if not os.path.isfile(prefix_path + ".reads"):
            print("No file", prefix_path + ".reads")
            return None
        if not os.path.isfile(prefix_path + ".ridx"):
            print("No file", prefix_path + ".ridx")
            return None
    return load_aindex(settings, skip_reads=skip_aindex, skip_aindex=skip_aindex)


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


def get_srandness(kmer, kmer2tf, k=23):
    ''' Wrapper that return number of + strand and - srand.
    '''
    poses = kmer2tf.pos(kmer)
    plus = 0
    minus = 0
    for pos in poses:
        _kmer = kmer2tf.reads[pos:pos+k]
        if kmer == _kmer:
            plus += 1
        else:
            minus += 1
    return plus, minus, len(poses)


def iter_reads_by_kmer(kmer, kmer2tf, used_reads=None, only_left=False, skip_multiple=True, k=23):
    ''' Yield 
        (start, next_read_start, read, pos_if_uniq|None, all_poses)

    - only_left: return only left reads
    - skip_multiple: skip if more then one hit in read

    '''

    rid2poses = kmer2tf.get_rid2poses(kmer)

    for rid in rid2poses:
        if used_reads and rid in used_reads:
            continue
        poses = rid2poses[rid]
        if skip_multiple:
            if len(poses) > 1:
                continue

        end = rid
        while True:
            if kmer2tf.end_cheker(end):
                break
            end += 1
        read = kmer2tf.reads[rid:end].decode("utf8")

        pos = poses[0]
        is_multiple_hit = len(poses) > 1
        print(read[pos:pos+k], kmer)
        if read[pos:pos+k] != kmer:
            read = get_revcomp(read)
            poses = [len(read) - x - k for x in poses]
            ori_pos = pos
            pos = poses[0]
            assert read[pos:pos+k] == kmer.encode("utf-8")
                
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
    hits = kmer2tf.get_rid2poses(kmer)

    for hit in hits:
        end = hit
        while True:
            if kmer2tf.end_cheker(end):
                break
            end += 1
        poses = hits[hit]
        read = kmer2tf.reads[hit:end].decode("utf-8")
        was_reversed = 0

        pos = poses[0]
        if read[pos:pos+k] != kmer:
            read = get_revcomp(read)
            poses = [len(read) - x - k for x in poses]
            pos = poses[0]
            was_reversed = 1
            print(read[pos:pos+k], kmer)
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


def get_left_right_distances(left_kmer, right_kmer, kmer2tf, k=23):
    '''
    Return a list of (rid, left_kmer_pos, right_kmer_pos) tuples.

    TODO: Handle cases: 1. distance in one subread; 2. distance in pair.
    TODO: implement it more efficently.
    '''
    hits = defaultdict(list)
    for pos in kmer2tf.pos(left_kmer):
        if kmer2tf.reads[pos:pos+k] == left_kmer.encode("utf-8"):
            start = kmer2tf.get_rid(pos)
            hits[start].append((0,pos))

    for pos in kmer2tf.pos(right_kmer):
        if kmer2tf.reads[pos:pos+k] == right_kmer.encode("utf-8"):
            start = kmer2tf.get_rid(pos)
            hits[start].append((1,pos))

    results = []
    for rid, hit in hits.items():
        if len(hit) == 1:
            continue
        both_kmers_found = set([x[0] for x in hit])
        if len(both_kmers_found) == 1:
            continue
        if len(hit) == 2:
            left_hit = hit[0][1]
            right_hit = hit[0][1]
            frag = kmer2tf.reads[min(left_hit, right_hit):max(left_hit, right_hit):]

            if b"~" in frag:
                results.append((rid, left_hit, right_hit, None))
            else:
                results.append((rid, left_hit, right_hit, right_hit-left_hit))
        else:
            results.append((rid, None, None, None))
        # for left_pos in [x[1] for x in hit if x[0] == 0]:
        #     for right_pos in [x[1] for x in hit if x[0] == 1]:
        #         results.append((rid, left_pos, right_pos))
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
        left, right = read.split("~")
        if pos < spring_pos:
            lefts.append("")
            rights.append(right)
            read = left
        else:
            lefts.append(left)
            rights.append("")
            pos = pos - len(left) - 1
            read = right
        max_pos = max(max_pos,pos)
        reads.append(read)
        starts.append(pos)
        rids.append(rid)
    max_length = max([len(x)+max_pos-starts[i] for i,x in enumerate(reads)])
    for i,read in enumerate(reads):
        separator = "N"
        reads[i] = separator*(max_pos-starts[i]) + read + separator * (max_length-max_pos+starts[i]-len(read))
    return max_pos, reads, lefts, rights, rids, starts

### Assembly-by-extension

def get_reads_for_assemby_by_kmer(kmer2tf, kmer, used_reads, compute_cov=True, k=23, mode=None):
    ''' Get reads prepared for assembly-by-extension. 
        Return sorted by pos list of (pos, read, rid, poses, cov)
        Mode: left, right
    '''    
    to_assembly = []
    for rid, poses in kmer2tf.get_rid2poses(kmer).items():
        if rid in used_reads:
            continue
        used_reads.add(rid)
        read = kmer2tf.get_read_by_rid(rid)

        spring_pos = None
        if mode:
            spring_pos = read.find("~")

        ori_poses = poses
        if not read[poses[0]:poses[0]+k] == kmer:
            read = get_revcomp(read)
            poses = [x for x in map(lambda x: len(read)-x-k, poses)][::-1]

        if mode == "left":
            read = read.split("~")[0]
            poses = [x for x in poses if x < spring_pos]
        elif mode == "right":
            read = read.split("~")[-1]
            poses = [x for x in poses if x > spring_pos]

        if not poses:
            continue

        cov = None
        if compute_cov:
            cov = [kmer2tf[read[i:i+k]] for i in range(len(read)-k+1)]
        to_assembly.append((poses[0], read, rid, ori_poses, cov))
    to_assembly.sort(reverse=True)
    return to_assembly