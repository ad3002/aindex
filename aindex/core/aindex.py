#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @created: 07.03.2015
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import os
import ctypes
from ctypes import cdll, c_void_p, c_char_p, c_uint64, c_uint32
import mmap
from collections import defaultdict
import importlib.resources as pkg_resources
from enum import IntEnum
from intervaltree import IntervalTree
from editdistance import eval as edit_distance
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# Note: Logging handlers should be configured in the main application.

# Load the shared library
with pkg_resources.path('aindex.core', 'python_wrapper.so') as dll_path:
    dll_path = str(dll_path)

if not os.path.exists(dll_path):
    logger.error(f"aindex's DLL was not found: {dll_path}")
    raise FileNotFoundError(f"aindex's DLL was not found: {dll_path}")

lib = cdll.LoadLibrary(dll_path)

# Define argument and return types for the shared library functions
lib.AindexWrapper_new.argtypes = []
lib.AindexWrapper_new.restype = c_void_p

lib.AindexWrapper_load.argtypes = [c_void_p, c_char_p, c_char_p]
lib.AindexWrapper_load.restype = None

lib.AindexWrapper_get.argtypes = [c_void_p, c_char_p]
lib.AindexWrapper_get.restype = c_uint64

lib.AindexWrapper_get_kid_by_kmer.argtypes = [c_void_p, c_char_p]
lib.AindexWrapper_get_kid_by_kmer.restype = c_uint64

lib.AindexWrapper_get_kmer_by_kid.argtypes = [c_void_p, c_uint64, c_char_p]
lib.AindexWrapper_get_kmer_by_kid.restype = None

lib.AindexWrapper_load_index.argtypes = [c_void_p, c_char_p, c_uint32]
lib.AindexWrapper_load_index.restype = None

lib.AindexWrapper_load_reads.argtypes = [c_void_p, c_char_p]
lib.AindexWrapper_load_reads.restype = None

lib.AindexWrapper_load_reads_index.argtypes = [c_void_p, c_char_p]
lib.AindexWrapper_load_reads_index.restype = None

lib.AindexWrapper_get_hash_size.argtypes = [c_void_p]
lib.AindexWrapper_get_hash_size.restype = c_uint64

lib.AindexWrapper_get_reads_size.argtypes = [c_void_p]
lib.AindexWrapper_get_reads_size.restype = c_uint64

lib.AindexWrapper_get_read.argtypes = [c_void_p, c_uint64, c_uint64, c_uint32]
lib.AindexWrapper_get_read.restype = c_char_p

lib.AindexWrapper_get_read_by_rid.argtypes = [c_void_p, c_uint64]
lib.AindexWrapper_get_read_by_rid.restype = c_char_p

lib.AindexWrapper_get_rid.argtypes = [c_void_p, c_uint64]
lib.AindexWrapper_get_rid.restype = c_uint64

lib.AindexWrapper_get_start.argtypes = [c_void_p, c_uint64]
lib.AindexWrapper_get_start.restype = c_uint64

lib.AindexWrapper_get_strand.argtypes = [c_void_p, c_char_p]
lib.AindexWrapper_get_strand.restype = c_uint64

lib.AindexWrapper_get_kmer.argtypes = [c_void_p, c_uint64, c_char_p, c_char_p]
lib.AindexWrapper_get_kmer.restype = c_uint64

lib.AindexWrapper_get_positions.argtypes = [c_void_p, ctypes.POINTER(c_uint64), c_char_p]
lib.AindexWrapper_get_positions.restype = None

lib.AindexWrapper_set_positions.argtypes = [c_void_p, ctypes.POINTER(c_uint64), c_char_p]
lib.AindexWrapper_set_positions.restype = None

class Strand(IntEnum):
    NOT_FOUND = 0
    FORWARD = 1
    REVERSE = 2

def get_revcomp(sequence: str) -> str:
    '''Return reverse complementary sequence.

    >>> get_revcomp('ATCGN')
    'NCGAT'

    '''
    complement = str.maketrans('ATCGNatcgn~[]', 'TAGCNtagcn~][')
    return sequence.translate(complement)[::-1]

def hamming_distance(s1: str, s2: str) -> int:
    """Get Hamming distance between two strings, ignoring positions with 'N'."""
    return sum(i != j for i, j in zip(s1, s2) if i != 'N' and j != 'N')

class AIndex:
    '''Wrapper for C++ AIndex implementation.'''

    def __init__(self, index_prefix: str):
        '''Initialize AIndex wrapper and load perfect hash.'''
        self.obj = lib.AindexWrapper_new()
        self.loaded_header = False
        self.loaded_intervals = False
        self.loaded_reads = False
        self.reads_size = 0

        required_files = [index_prefix + ext for ext in [".pf", ".tf.bin", ".kmers.bin"]]
        missing_files = [f for f in required_files if not os.path.isfile(f)]
        if missing_files:
            logger.error(f"One or more index files were not found: {', '.join(missing_files)}")
            raise FileNotFoundError(f"One or more index files were not found: {', '.join(missing_files)}")

        tf_file = index_prefix + ".tf.bin"
        lib.AindexWrapper_load(self.obj, index_prefix.encode('utf-8'), tf_file.encode('utf-8'))

    def load(self, index_prefix: str, max_tf: int):
        '''Load AIndex with a maximum term frequency limit.'''
        logger.info(f"Loading AIndex: {index_prefix}.*")

        required_files = [index_prefix + ext for ext in [".pf", ".tf.bin", ".kmers.bin", ".index.bin", ".indices.bin", ".pos.bin"]]
        missing_files = [f for f in required_files if not os.path.isfile(f)]
        if missing_files:
            logger.error(f"One or more index files were not found: {', '.join(missing_files)}")
            raise FileNotFoundError(f"One or more index files were not found: {', '.join(missing_files)}")

        self.max_tf = max_tf

        tf_file = index_prefix + ".tf.bin"
        lib.AindexWrapper_load_index(self.obj, index_prefix.encode('utf-8'), c_uint32(max_tf))

    def load_reads(self, reads_file: str):
        '''Load reads using the shared library.'''
        if not os.path.isfile(reads_file):
            logger.error(f"Reads file was not found: {reads_file}")
            raise FileNotFoundError(f"Reads file was not found: {reads_file}")

        logger.info(f"Loading reads: {reads_file}")
        lib.AindexWrapper_load_reads(self.obj, reads_file.encode('utf-8'))
        self.reads_size = lib.AindexWrapper_get_reads_size(self.obj)
        logger.info(f"Loaded {self.reads_size} reads.")

    def load_reads_index(self, index_file: str, header_file: str = None):
        '''Load reads index and optional headers.'''
        logger.info(f"Loading reads index: {index_file}")
        self.rid2start = {}
        self.IT = IntervalTree()
        self.chrm2start = {}
        self.headers = {}

        with open(index_file) as fh:
            for line in fh:
                rid_str, start_str, end_str = line.strip().split("\t")
                rid = int(rid_str)
                start = int(start_str)
                end = int(end_str)
                self.rid2start[rid] = (start, end)
                self.IT.addi(start, end, rid)
        self.loaded_intervals = True

        if header_file:
            logger.info(f"Loading headers: {header_file}")
            with open(header_file) as fh:
                for rid, line in enumerate(fh):
                    head, start_str, length_str = line.strip().split("\t")
                    start = int(start_str)
                    length = int(length_str)
                    self.headers[rid] = head
                    chrm = head.split()[0].split(".")[0]
                    self.chrm2start[chrm] = start
                    self.IT.addi(start, start + length, head)
            self.loaded_header = True

    def get_hash_size(self) -> int:
        '''Get hash size.'''
        return lib.AindexWrapper_get_hash_size(self.obj)

    def __len__(self):
        '''Get number of kmers.'''
        return self.get_hash_size()

    def __getitem__(self, kmer: str) -> int:
        '''Return term frequency for kmer.'''
        return lib.AindexWrapper_get(self.obj, kmer.encode('utf-8'))

    def get_strand(self, kmer: str) -> Strand:
        '''Return strand for kmer.'''
        result = lib.AindexWrapper_get_strand(self.obj, kmer.encode('utf-8'))
        return Strand(result)

    def get_kid_by_kmer(self, kmer: str) -> int:
        '''Return kmer ID for kmer.'''
        return lib.AindexWrapper_get_kid_by_kmer(self.obj, kmer.encode('utf-8'))

    def get_kmer_by_kid(self, kid: int, k: int = 23) -> str:
        '''Return kmer by kmer ID.'''
        kmer = ctypes.create_string_buffer(k)
        lib.AindexWrapper_get_kmer_by_kid(self.obj, c_uint64(kid), kmer)
        return kmer.value.decode('utf-8')

    def get_kmer_info_by_kid(self, kid: int, k: int = 23):
        '''Get kmer, reverse complement kmer, and corresponding term frequency for a given kmer ID.'''
        kmer = ctypes.create_string_buffer(k)
        rkmer = ctypes.create_string_buffer(k)
        tf = lib.AindexWrapper_get_kmer(self.obj, c_uint64(kid), kmer, rkmer)
        return kmer.value.decode('utf-8'), rkmer.value.decode('utf-8'), tf

    def __contains__(self, kmer: str) -> bool:
        """Check if a kmer exists in the index."""
        return self[kmer] > 0

    def get(self, kmer: str, default: int = 0) -> int:
        """Get term frequency for kmer, return default if not found."""
        tf = self[kmer]
        return tf if tf > 0 else default

    def get_rid(self, pos: int) -> int:
        '''Get read ID by position in read file.'''
        return lib.AindexWrapper_get_rid(self.obj, c_uint64(pos))

    def get_start(self, pos: int) -> int:
        '''Get start position of read by position in read file.'''
        return lib.AindexWrapper_get_start(self.obj, c_uint64(pos))

    def get_read_by_rid(self, rid: int) -> str:
        '''Get read sequence as string by read ID.'''
        read = lib.AindexWrapper_get_read_by_rid(self.obj, c_uint64(rid))
        return read.decode('utf-8') if read else ''

    def get_read(self, start: int, end: int, revcomp: bool = False) -> str:
        '''Get read by start and end positions.'''
        read = lib.AindexWrapper_get_read(self.obj, c_uint64(start), c_uint64(end), int(revcomp))
        return read.decode('utf-8') if read else ''

    def iter_reads(self):
        '''Iterate over reads and yield (read_id, read).'''
        if self.reads_size == 0:
            logger.error("Reads were not loaded.")
            raise RuntimeError("Reads were not loaded.")

        for rid in range(self.reads_size):
            yield rid, self.get_read_by_rid(rid)

    def iter_reads_se(self):
        '''Iterate over reads and yield (read_id, subread_index, subread).'''
        if self.reads_size == 0:
            logger.error("Reads were not loaded.")
            raise RuntimeError("Reads were not loaded.")

        for rid in range(self.reads_size):
            read = self.get_read_by_rid(rid)
            subreads = read.split("~")
            for idx, subread in enumerate(subreads):
                yield rid, idx, subread

    def pos(self, kmer: str) -> list:
        '''Return list of positions for a given kmer.'''
        n = self.max_tf
        r = (c_uint64 * n)()
        lib.AindexWrapper_get_positions(self.obj, r, kmer.encode('utf-8'))
        return [r[i] - 1 for i in range(n) if r[i] > 0]

    def get_header(self, pos: int) -> str:
        '''Get header information for a position.'''
        if not self.loaded_header:
            return None
        intervals = self.IT[pos]
        if intervals:
            rid = next(iter(intervals)).data
            return self.headers.get(rid, '')
        return ''

    def iter_sequence_kmers(self, sequence: str, k: int = 23):
        '''Iterate over kmers in a sequence and yield (kmer, term_frequency).'''
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            if '\n' in kmer or '~' in kmer:
                continue
            yield kmer, self[kmer]

    def get_sequence_coverage(self, seq: str, cutoff: int = 0, k: int = 23) -> list:
        '''Get coverage of a sequence based on kmers.'''
        coverage = [0] * len(seq)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            tf = self[kmer]
            if tf >= cutoff:
                coverage[i] = tf
        return coverage

    def print_sequence_coverage(self, seq: str, cutoff: int = 0):
        '''Print sequence coverage and return list of term frequencies for each kmer.'''
        coverage = self.get_sequence_coverage(seq, cutoff)
        for i, tf in enumerate(coverage):
            kmer = seq[i:i + 23]
            print(f"{i}\t{kmer}\t{tf}")
        return coverage

    def get_rid2poses(self, kmer: str) -> dict:
        '''Return a mapping from read ID to positions in read for a given kmer.'''
        poses = self.pos(kmer)
        hits = defaultdict(list)
        for pos in poses:
            rid = self.get_rid(pos)
            start = self.get_start(pos)
            hits[rid].append(pos - start)
        return hits

    def set(self, poses_array: list, kmer: str, batch_size: int):
        '''Update kmer batch (this function is not fully implemented).'''
        logger.warning("WARNING: called a function with the fixed batch size.")
        n = batch_size * 2
        r = (c_uint64 * n)()
        for i, (rid, pos) in enumerate(poses_array):
            r[i + batch_size] = c_uint64(rid)
            r[i] = c_uint64(pos)
        lib.AindexWrapper_set_positions(self.obj, r, kmer.encode('utf-8'))

def get_aindex(prefix_path: str, skip_aindex: bool = False, max_tf: int = 1_000_000) -> AIndex:
    required_files = [
        f"{prefix_path}.23.pf",
        f"{prefix_path}.23.tf.bin",
        f"{prefix_path}.23.kmers.bin",
    ]
    if not skip_aindex:
        required_files.extend([
            f"{prefix_path}.23.index.bin",
            f"{prefix_path}.23.indices.bin",
            f"{prefix_path}.23.pos.bin",
            f"{prefix_path}.reads",
            f"{prefix_path}.ridx",
        ])

    missing_files = [file for file in required_files if not os.path.isfile(file)]
    if missing_files:
        logger.error(f"Required files not found: {', '.join(missing_files)}")
        raise FileNotFoundError(f"Required files not found: {', '.join(missing_files)}")

    settings = {
        "index_prefix": f"{prefix_path}.23",
        "aindex_prefix": f"{prefix_path}.23",
        "reads_file": f"{prefix_path}.reads",
        "max_tf": max_tf,
    }

    return load_aindex(settings, skip_reads=skip_aindex, skip_aindex=skip_aindex)

def load_aindex(settings: dict, prefix: str = None, reads: str = None, aindex_prefix: str = None, skip_reads: bool = False, skip_aindex: bool = False) -> AIndex:
    '''Load AIndex, including optional reads and aindex data.'''
    if prefix is None:
        prefix = settings.get("index_prefix")
    if reads is None and not skip_reads:
        reads = settings.get("reads_file")

    max_tf = settings.get("max_tf", 10000)
    if aindex_prefix is None and not skip_aindex:
        aindex_prefix = settings.get("aindex_prefix")

    aindex = AIndex(prefix)
    aindex.max_tf = max_tf

    if not skip_reads:
        aindex.load_reads(reads)
    if not skip_aindex:
        aindex.load(aindex_prefix, max_tf)

    return aindex

def get_srandness(kmer: str, aindex: AIndex, k: int = 23) -> tuple:
    '''Return the number of plus and minus strand occurrences for a kmer.'''
    poses = aindex.pos(kmer)
    plus = minus = 0
    for pos in poses:
        read_kmer = aindex.get_read(pos, pos + k)
        if kmer == read_kmer:
            plus += 1
        else:
            minus += 1
    return plus, minus, len(poses)

def iter_reads_by_kmer(kmer: str, aindex: AIndex, used_reads: set = None, skip_multiple: bool = False, k: int = 23):
    '''Yield reads containing the given kmer.'''
    if used_reads is None:
        used_reads = set()
    rid2poses = aindex.get_rid2poses(kmer)
    for rid, poses in rid2poses.items():
        if rid in used_reads:
            continue
        used_reads.add(rid)
        read = aindex.get_read_by_rid(rid)
        if skip_multiple and len(poses) > 1:
            continue
        for pos in poses:
            if not read[pos:pos + k] == kmer:
                read = get_revcomp(read)
                poses = [len(read) - x - k for x in poses]
                pos = poses[0]
            yield [rid, pos, read, poses]

def iter_reads_by_sequence(sequence: str, aindex: AIndex, hd: int = None, ed: int = None, used_reads: set = None, skip_multiple: bool = False, k: int = 23):
    '''Yield reads containing the given sequence.'''
    if len(sequence) >= k:
        kmer = sequence[:k]
        n = len(sequence)
        for rid, pos, read, poses in iter_reads_by_kmer(kmer, aindex, used_reads=used_reads, skip_multiple=skip_multiple, k=k):
            for pos in poses:
                if not hd and sequence in read:
                    yield rid, pos, read, poses
                elif hd:
                    fragment = read[pos:pos + n]
                    if len(fragment) == n and hamming_distance(fragment, sequence) <= hd:
                        yield rid, pos, read, poses, hd
                elif ed:
                    fragment = read[pos:pos + n]
                    if len(fragment) == n and edit_distance(fragment, sequence) <= ed:
                        yield rid, pos, read, poses, ed
    else:
        yield None

def iter_reads_se_by_kmer(kmer: str, aindex: AIndex, used_reads: set = None, k: int = 23):
    '''Yield split reads containing the given kmer.'''
    if used_reads is None:
        used_reads = set()
    for rid, pos, read, poses in iter_reads_by_kmer(kmer, aindex, used_reads=used_reads, k=k):
        spring_pos = read.find("~")
        if spring_pos == -1:
            yield [rid, pos, read, -1]
            continue
        left, right = read.split("~")
        if pos < spring_pos:
            read = left
            pos = pos
            yield [rid, pos, read, 0]
        else:
            read = right
            pos = pos - spring_pos - 1
            yield [rid, pos, read, 1]

def get_left_right_distances(left_kmer: str, right_kmer: str, aindex: AIndex, k: int = 23):
    '''Return distances between left and right kmers in reads.'''
    hits = defaultdict(list)
    for pos in aindex.pos(left_kmer):
        rid = aindex.get_rid(pos)
        hits[rid].append((0, pos))
    for pos in aindex.pos(right_kmer):
        rid = aindex.get_rid(pos)
        hits[rid].append((1, pos))
    for rid, hit_list in hits.items():
        if len(hit_list) != 2:
            continue
        hit_list.sort()
        (type1, pos1), (type2, pos2) = hit_list
        if type1 == type2:
            continue
        start, end = sorted([pos1, pos2])
        fragment = aindex.get_read(start, end + k)
        reversed_read = False
        if not fragment:
            fragment = aindex.get_read(end, start + k, revcomp=True)
            reversed_read = True
        has_spring = '~' in fragment
        yield rid, start, end, len(fragment), fragment, has_spring, reversed_read

def get_layout_from_reads(kmer: str, aindex: AIndex, used_reads: set = None, k: int = 23, space: str = "N"):
    '''Get layout and flanking reads for a given kmer.'''
    if used_reads is None:
        used_reads = set()
    max_pos = 0
    reads = []
    lefts = []
    rights = []
    rids = []
    starts = []
    for rid, pos, read, poses in iter_reads_by_kmer(kmer, aindex, used_reads, k=k):
        spring_pos = read.find("~")
        if spring_pos > -1:
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
        max_pos = max(max_pos, pos)
        reads.append(read)
        starts.append(pos)
        rids.append(rid)
    max_length = max(len(read) + max_pos - starts[i] for i, read in enumerate(reads))
    aligned_reads = [space * (max_pos - starts[i]) + read + space * (max_length - (max_pos - starts[i]) - len(read)) for i, read in enumerate(reads)]
    return max_pos, aligned_reads, lefts, rights, rids, starts

### Assembly-by-extension

# def get_reads_for_assemby_by_kmer(kmer2tf, kmer, used_reads, compute_cov=True, k=23, mode=None):
#     ''' Get reads prepared for assembly-by-extension.
#         Return sorted by pos list of (pos, read, rid, poses, cov)
#         Mode: left, right
#     '''
#     to_assembly = []
#     for rid, poses in kmer2tf.get_rid2poses(kmer).items():
#         if rid in used_reads:
#             continue
#         used_reads.add(rid)
#         read = kmer2tf.get_read_by_rid(rid)

#         spring_pos = None
#         if mode:
#             spring_pos = read.find("~")

#         ori_poses = poses
#         if not read[poses[0]:poses[0]+k] == kmer:
#             read = get_revcomp(read)
#             poses = [x for x in map(lambda x: len(read)-x-k, poses)][::-1]

#         if mode == "left":
#             read = read.split("~")[0]
#             poses = [x for x in poses if x < spring_pos]
#         elif mode == "right":
#             read = read.split("~")[-1]
#             poses = [x for x in poses if x > spring_pos]

#         if not poses:
#             continue

#         cov = None
#         if compute_cov:
#             cov = [kmer2tf[read[i:i+k]] for i in range(len(read)-k+1)]
#         to_assembly.append((poses[0], read, rid, ori_poses, cov))
#     to_assembly.sort(reverse=True)
#     return to_assembly
