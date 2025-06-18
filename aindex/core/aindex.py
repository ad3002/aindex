#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @created: 07.03.2015
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

"""
Aindex Python API using pybind11 bindings
Modern, safe, and Pythonic interface for k-mer indexing
All file paths must be provided explicitly - no automatic name generation
"""

from typing import Optional, Dict, Any, List, Tuple, Iterator
try:
    from . import aindex_cpp
except ImportError:
    import aindex_cpp
import os
from collections import defaultdict
from enum import IntEnum
from intervaltree import IntervalTree
from editdistance import eval as edit_distance
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class Strand(IntEnum):
    NOT_FOUND = 0
    FORWARD = 1
    REVERSE = 2

def get_revcomp(sequence: str) -> str:
    """Return reverse complementary sequence.

    >>> get_revcomp('ATCGN')
    'NCGAT'

    """
    complement = str.maketrans('ATCGNatcgn~[]', 'TAGCNtagcn~][')
    return sequence.translate(complement)[::-1]

def hamming_distance(s1: str, s2: str) -> int:
    """Get Hamming distance between two strings, ignoring positions with 'N'."""
    return sum(i != j for i, j in zip(s1, s2) if i != 'N' and j != 'N')

class AIndex:
    """
    Python wrapper for the Aindex C++ library using pybind11
    All file paths must be provided explicitly - no automatic name generation
    """
    
    def __init__(self):
        self._wrapper = aindex_cpp.AindexWrapper()
        self._loaded = False
        self.reads_size = 0
        self.loaded_header = False
        self.loaded_intervals = False
        self.loaded_reads = False
        self.max_tf = 0
        
    def load_hash(self, hash_file: str, tf_file: str, kmers_bin_file: str, kmers_text_file: str):
        """
        Load hash index from explicit file paths
        All files must exist and be provided explicitly
        """
        # Validate required files exist
        for fname, fpath in [('hash', hash_file), ('tf', tf_file), 
                           ('kmers_bin', kmers_bin_file)]:
            if not os.path.exists(fpath):
                raise FileNotFoundError(f"{fname} file not found: {fpath}")
        
        # Only check kmers_text_file if it's not empty
        if kmers_text_file and not os.path.exists(kmers_text_file):
            raise FileNotFoundError(f"kmers_text file not found: {kmers_text_file}")
        
        self._wrapper.load(hash_file, tf_file, kmers_bin_file, kmers_text_file)
        self._loaded = True
    
    def load_hash_file(self, hash_file: str, tf_file: str, kmers_bin_file: str, kmers_text_file: str):
        """
        Load hash file with explicit paths (alias for load_hash)
        """
        self.load_hash(hash_file, tf_file, kmers_bin_file, kmers_text_file)
        
    def load_reads(self, reads_file: str):
        """Load reads file"""
        if not os.path.exists(reads_file):
            raise FileNotFoundError(f"Reads file not found: {reads_file}")
        self._wrapper.load_reads(reads_file)
        self.reads_size = self._wrapper.reads_size
        
    def load_aindex(self, index_file: str, indices_file: str, max_tf: int):
        """Load aindex files with explicit paths"""
        for fname, fpath in [('index', index_file), ('indices', indices_file)]:
            if not os.path.exists(fpath):
                raise FileNotFoundError(f"{fname} file not found: {fpath}")
        self._wrapper.load_aindex(index_file, indices_file, max_tf)
        
    def load_reads_index(self, index_file: str, header_file: str = None):
        """Load reads index and optional headers."""
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

    def get_tf_value(self, kmer: str) -> int:
        """Get term frequency for a kmer"""
        if not self._loaded:
            return 0  # Return 0 if no index is loaded
        return self._wrapper.get_tf_value(kmer)
        
    def get_tf_values(self, kmers: List[str]) -> List[int]:
        """Get term frequencies for multiple kmers"""
        if not self._loaded:
            return [0] * len(kmers)  # Return zeros if no index is loaded
        return self._wrapper.get_tf_values(kmers)
        
    def get_tf_values_13mer(self, kmers: List[str]) -> List[int]:
        """Get term frequency values for 13-mers using batch processing"""
        if not self._loaded:
            return [0] * len(kmers)
        return self._wrapper.get_tf_values_13mer(kmers)
        
    def get_hash_value(self, kmer: str) -> int:
        """Get hash value for a kmer"""
        if not self._loaded:
            raise RuntimeError("Index not loaded")
        return self._wrapper.get_hash_value(kmer)
        
    def get_hash_values(self, kmers: List[str]) -> List[int]:
        """Get hash values for multiple kmers"""
        if not self._loaded:
            raise RuntimeError("Index not loaded")
        return self._wrapper.get_hash_values(kmers)
        
    def get_reads_by_kmer(self, kmer: str, max_reads: int = 100) -> List[str]:
        """Get reads containing a specific kmer"""
        if not self._wrapper.aindex_loaded:
            raise RuntimeError("Aindex not loaded")
        return self._wrapper.get_reads_se_by_kmer(kmer, max_reads)
        
    def get_read_by_rid(self, rid: int) -> str:
        """Get read by read ID"""
        return self._wrapper.get_read_by_rid(rid)
        
    def get_read(self, start: int, end: int, revcomp: bool = False) -> str:
        """Get read by start and end positions"""
        return self._wrapper.get_read(start, end, revcomp)
        
    def get_kid_by_kmer(self, kmer: str) -> int:
        """Get kmer ID by kmer"""
        if not self._loaded:
            raise RuntimeError("Index not loaded")
        return self._wrapper.get_kid_by_kmer(kmer)
        
    def get_kmer_by_kid(self, kid: int) -> str:
        """Get kmer by kmer ID"""
        if not self._loaded:
            raise RuntimeError("Index not loaded")
        return self._wrapper.get_kmer_by_kid(kid)
        
    def get_strand(self, kmer: str) -> Strand:
        """Return strand for kmer."""
        if not self._loaded:
            raise RuntimeError("Index not loaded")
        result = self._wrapper.get_strand(kmer)
        return Strand(result)
        
    def get_kmer_info(self, kid: int) -> Tuple[str, str, int]:
        """Get kmer info by kmer ID (kmer, rkmer, tf)"""
        if not self._loaded:
            raise RuntimeError("Index not loaded")
            
        # C++/Python string references don't work as expected for output parameters
        # We need to get the kmer directly from the kid using get_kmer_by_kid
        # and then compute its reverse complement for the return values
        kmer = self.get_kmer_by_kid(kid)
        from Bio.Seq import Seq
        rkmer = str(Seq(kmer).reverse_complement())
        tf = self.get_tf_value(kmer)
        return kmer, rkmer, tf
        
    def get_rid(self, pos: int) -> int:
        """Get read ID by position"""
        if not self._wrapper.aindex_loaded:
            raise RuntimeError("Aindex not loaded")
        return self._wrapper.get_rid(pos)
        
    def get_start(self, pos: int) -> int:
        """Get start position by position"""
        if not self._wrapper.aindex_loaded:
            raise RuntimeError("Aindex not loaded")
        return self._wrapper.get_start(pos)
        
    def get_positions(self, kmer: str) -> List[int]:
        """Get positions for kmer (supports both 13-mers and 23-mers)"""
        if len(kmer) == 13:
            # For 13-mers, require 13-mer aindex to be loaded
            # The C++ wrapper will handle routing to the appropriate function
            return self._wrapper.get_positions(kmer)
        elif len(kmer) == 23:
            # For 23-mers, require traditional aindex to be loaded
            if not self._wrapper.aindex_loaded:
                raise RuntimeError("23-mer Aindex not loaded")
            return self._wrapper.get_positions(kmer)
        else:
            raise ValueError(f"Unsupported k-mer length: {len(kmer)}. Only 13-mers and 23-mers are supported.")
        
    def get_positions_13mer(self, kmer: str) -> List[int]:
        """Get positions for 13-mers (direct method)"""
        return self._wrapper.get_positions_13mer(kmer)
        
    def get_hash_size(self) -> int:
        """Get hash size"""
        if not self._loaded:
            raise RuntimeError("Index not loaded")
        return self._wrapper.get_hash_size()
        
    def get_reads_size(self) -> int:
        """Get reads size"""
        return self._wrapper.get_reads_size()
        
    def __len__(self) -> int:
        """Get number of kmers"""
        return self.get_hash_size()
        
    def __getitem__(self, kmer: str) -> int:
        """Return term frequency for kmer"""
        return self.get_tf_value(kmer)
        
    def __contains__(self, kmer: str) -> bool:
        """Check if a kmer exists in the index"""
        return self[kmer] > 0
        
    def get(self, kmer: str, default: int = 0) -> int:
        """Get term frequency for kmer, return default if not found"""
        tf = self[kmer]
        return tf if tf > 0 else default

    def get_kmer_info_by_kid(self, kid: int, k: int = 23):
        """Get kmer, reverse complement kmer, and corresponding term frequency for a given kmer ID."""
        kmer, rkmer, tf = self.get_kmer_info(kid)
        return kmer, rkmer, tf

    def iter_reads(self):
        """Iterate over reads and yield (read_id, read)."""
        if self.reads_size == 0:
            logger.error("Reads were not loaded.")
            raise RuntimeError("Reads were not loaded.")

        for rid in range(self.n_reads):
            yield rid, self.get_read_by_rid(rid)

    def iter_reads_se(self):
        """Iterate over reads and yield (read_id, subread_index, subread)."""
        if self.reads_size == 0:
            logger.error("Reads were not loaded.")
            raise RuntimeError("Reads were not loaded.")

        for rid in range(self.n_reads):
            read = self.get_read_by_rid(rid)
            subreads = read.split("~")
            for idx, subread in enumerate(subreads):
                yield rid, idx, subread

    def pos(self, kmer: str) -> list:
        """Return list of positions for a given kmer."""
        return self.get_positions(kmer)

    def get_header(self, pos: int) -> str:
        """Get header information for a position."""
        if not self.loaded_header:
            return None
        intervals = self.IT[pos]
        if intervals:
            rid = next(iter(intervals)).data
            return self.headers.get(rid, '')
        return ''

    def iter_sequence_kmers(self, sequence: str, k: int = 23):
        """Iterate over kmers in a sequence and yield (kmer, term_frequency)."""
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            if '\n' in kmer or '~' in kmer:
                continue
            yield kmer, self[kmer]

    def get_sequence_coverage(self, seq: str, cutoff: int = 0, k: int = 23) -> list:
        """Get coverage of a sequence based on kmers."""
        coverage = [0] * (len(seq) - k + 1)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            tf = self[kmer]
            if tf >= cutoff:
                coverage[i] = tf
        return coverage

    def print_sequence_coverage(self, seq: str, cutoff: int = 0):
        """Print sequence coverage and return list of term frequencies for each kmer."""
        coverage = self.get_sequence_coverage(seq, cutoff)
        for i, tf in enumerate(coverage):
            # Only print complete k-mers (avoid printing shorter k-mers at the end)
            kmer = seq[i:i + 23]
            print(f"{i}\t{kmer}\t{tf}")
        return coverage

    def get_rid2poses(self, kmer: str) -> dict:
        """Return a mapping from read ID to positions in read for a given kmer."""
        poses = self.pos(kmer)
        hits = defaultdict(list)
        for pos in poses:
            rid = self.get_rid(pos)
            start = self.get_start(pos)
            hits[rid].append(pos - start)
        return hits

        
    @property
    def n_reads(self) -> int:
        """Number of reads"""
        return self._wrapper.n_reads
        
    @property 
    def n_kmers(self) -> int:
        """Number of kmers"""
        return self._wrapper.n_kmers
        
    @property
    def aindex_loaded(self) -> bool:
        """Whether aindex is loaded"""
        return self._wrapper.aindex_loaded

    def load_13mer_index(self, hash_file: str, tf_file: str):
        """
        Load 13-mer index (hash and term frequencies) from explicit file paths
        """
        if not os.path.exists(hash_file):
            raise FileNotFoundError(f"13-mer hash file not found: {hash_file}")
        if not os.path.exists(tf_file):
            raise FileNotFoundError(f"13-mer tf file not found: {tf_file}")
        
        self._wrapper.load_13mer_index(hash_file, tf_file)
        self._loaded = True
        
    def load_13mer_aindex(self, index_file: str, indices_file: str):
        """
        Load 13-mer position index from explicit file paths
        """
        for fname, fpath in [('index', index_file), ('indices', indices_file)]:
            if not os.path.exists(fpath):
                raise FileNotFoundError(f"13-mer {fname} file not found: {fpath}")
        
        self._wrapper.load_13mer_aindex(index_file, indices_file)
        
    @staticmethod
    def load_13mer_index_static(hash_file: str, tf_file: str) -> 'AIndex':
        """
        Load 13-mer index with automatic mode detection
    
        Args:
            hash_file: Path to .pf file (from build_13mer_hash)
            tf_file: Path to .tf.bin file (from count_kmers13)
    
        Returns:
            AIndex instance configured for 13-mer mode
            
        Example:
            >>> index = load_13mer_index('13mers.pf', '13mers.tf.bin')
            >>> tf_value = index.get_tf_value('ATCGATCGATCGA')
            >>> all_counts = index.get_13mer_tf_array()
        """
        index = AIndex()
        index.load_13mer_index(hash_file, tf_file)
        return index

    @staticmethod
    def load_23mer_index(hash_file: str, tf_file: str, kmers_bin_file: str, kmers_text_file: str = "") -> 'AIndex':
        """
        Load traditional 23-mer index
    
        Args:
            hash_file: Path to hash file
            tf_file: Path to tf file  
            kmers_bin_file: Path to binary kmers file
            kmers_text_file: Path to text kmers file (optional)
    
        Returns:
            AIndex instance configured for 23-mer mode
        """
        index = AIndex()
        index.load_hash(hash_file, tf_file, kmers_bin_file, kmers_text_file)
        return index

    @staticmethod
    def load_from_prefix(prefix: str, kmer_size: Optional[int] = None, max_tf: int = 100000, 
                        load_aindex: bool = True, load_reads: bool = False) -> 'AIndex':
        """
        Load index from prefix with auto-detection or manual specification
        
        Args:
            prefix: File prefix (e.g., 'reads.23', 'reads.13', 'mydata')
            kmer_size: Manual k-mer size specification (13 or 23). If None, auto-detect from files
            max_tf: Maximum term frequency for aindex loading (23-mers only)
            load_aindex: Whether to load position index (default: True)
            load_reads: Whether to load reads file (default: False)
            
        Returns:
            AIndex instance configured for detected/specified mode
            
        Examples:
            >>> # Auto-detect from files
            >>> index = AIndex.load_from_prefix('reads.23')
            >>> index = AIndex.load_from_prefix('reads.13')
            
            >>> # Manual specification
            >>> index = AIndex.load_from_prefix('mydata', kmer_size=23, max_tf=100)
            >>> index = AIndex.load_from_prefix('mydata', kmer_size=13)
            
            >>> # Load without aindex
            >>> index = AIndex.load_from_prefix('reads.23', load_aindex=False)
            
            >>> # Load with reads
            >>> index = AIndex.load_from_prefix('reads.23', load_reads=True)
        """
        index = AIndex()
        
        # Auto-detect k-mer size if not specified
        if kmer_size is None:
            # Check for 13-mer files
            hash_13_file = f"{prefix}.pf"
            tf_13_file = f"{prefix}.tf.bin"
            
            # Check for 23-mer files  
            pf_23_file = f"{prefix}.pf"
            tf_23_file = f"{prefix}.tf.bin"
            kmers_bin_file = f"{prefix}.kmers.bin"
            
            if os.path.exists(hash_13_file) and os.path.exists(tf_13_file):
                kmer_size = 13
                logger.info(f"Auto-detected 13-mer mode from files: {hash_13_file}, {tf_13_file}")
            elif os.path.exists(pf_23_file) and os.path.exists(tf_23_file) and os.path.exists(kmers_bin_file):
                kmer_size = 23
                logger.info(f"Auto-detected 23-mer mode from files: {pf_23_file}, {tf_23_file}, {kmers_bin_file}")
            else:
                raise FileNotFoundError(
                    f"Could not auto-detect k-mer size for prefix '{prefix}'. "
                    f"Expected either:\n"
                    f"  13-mer: {hash_13_file} + {tf_13_file}\n"
                    f"  23-mer: {pf_23_file} + {tf_23_file} + {kmers_bin_file}"
                )
        
        # Determine reads file path if loading reads
        reads_file = ""
        if load_reads:
            reads_file = f"{prefix}.reads"
            if not os.path.exists(reads_file):
                reads_file = reads_file.replace(".23.", ".")
                reads_file = reads_file.replace(".13.", ".")
                if not os.path.exists(reads_file):
                    logger.warning(f"Reads file not found: {reads_file}")
                    reads_file = ""
        
        # Load based on kmer size
        if kmer_size == 13:
            index.load_from_prefix_13mer(prefix, load_aindex=load_aindex, reads_file=reads_file)
        elif kmer_size == 23:
            if max_tf is None:
                max_tf = 100000  # Default value for 23-mers
            index.load_from_prefix_23mer(prefix, max_tf=max_tf, load_aindex=load_aindex, reads_file=reads_file)
        else:
            raise ValueError(f"Unsupported kmer size: {kmer_size}. Only 13 and 23 are supported.")
        
        return index

    def load_from_prefix_23mer(self, prefix: str, max_tf: int = 100, load_aindex: bool = True, reads_file: str = ""):
        """
        Load 23-mer index and optionally aindex from prefix
        
        Args:
            prefix: File prefix
            max_tf: Maximum term frequency for aindex loading
            load_aindex: Whether to load position index
            reads_file: Optional path to reads file (if not provided, no reads will be loaded)
        """
        # Load hash index
        self._wrapper.load_from_prefix_23mer(prefix, reads_file)
        self._loaded = True
        
        # Load aindex if requested
        if load_aindex:
            try:
                self._wrapper.load_aindex_from_prefix_23mer(prefix, max_tf, reads_file)
                logger.info(f"23-mer AIndex loaded from prefix: {prefix}")
            except Exception as e:
                logger.warning(f"Could not load 23-mer AIndex from prefix {prefix}: {e}")

    def load_from_prefix_13mer(self, prefix: str, load_aindex: bool = True, reads_file: str = ""):
        """
        Load 13-mer index and optionally aindex from prefix
        
        Args:
            prefix: File prefix
            load_aindex: Whether to load position index
            reads_file: Optional path to reads file (if not provided, no reads will be loaded)
        """
        # Load hash index
        self._wrapper.load_from_prefix_13mer(prefix, reads_file)
        self._loaded = True
        
        # Load aindex if requested
        if load_aindex:
            try:
                self._wrapper.load_aindex_from_prefix_13mer(prefix, reads_file)
                logger.info(f"13-mer AIndex loaded from prefix: {prefix}")
            except Exception as e:
                logger.warning(f"Could not load 13-mer AIndex from prefix {prefix}: {e}")

    def get_13mer_tf_array(self) -> List[int]:
        """
        Get direct access to 13-mer tf array
        
        Returns:
            List of term frequencies for all possible 13-mers (4^13 = 67,108,864 elements)
        """
        return self._wrapper.get_13mer_tf_array()

    def get_tf_by_index_13mer(self, index: int) -> int:
        """
        Get tf value by direct array index for 13-mers
        
        Args:
            index: Array index (0 to 4^13-1)
            
        Returns:
            Term frequency value
        """
        return self._wrapper.get_tf_by_index_13mer(index)

    def get_index_info(self) -> str:
        """
        Get statistics about loaded index
        
        Returns:
            String with index information including mode, k-mer counts, etc.
        """
        return self._wrapper.get_index_info()

    def _index_to_13mer(self, index: int) -> str:
        """
        Convert index to 13-mer string
        
        Args:
            index: Index in the range [0, 4^13-1]
            
        Returns:
            13-mer string
        """
        nucleotides = ['A', 'C', 'G', 'T']
        kmer = []
        temp_index = index
        
        for i in range(13):
            kmer.append(nucleotides[temp_index % 4])
            temp_index //= 4
            
        return ''.join(reversed(kmer))

    def iter_kmers_by_frequency(self, min_tf: int = 1, max_kmers: Optional[int] = None, 
                               kmer_type: str = "auto") -> Iterator[Tuple[str, int]]:
        """
        Iterate over k-mers sorted by term frequency (most frequent first)
        
        Args:
            min_tf: Minimum term frequency threshold (default: 1)
            max_kmers: Maximum number of k-mers to return (default: None - return all)
            kmer_type: Type of k-mers to iterate over: "13mer", "23mer", or "auto" (default: "auto")
            
        Yields:
            Tuple of (kmer, term_frequency) sorted by decreasing frequency
            
        Examples:
            >>> # Get top 100 most frequent 13-mers
            >>> for kmer, tf in index.iter_kmers_by_frequency(max_kmers=100, kmer_type="13mer"):
            ...     print(f"{kmer}: {tf}")
            
            >>> # Get all 13-mers with tf >= 10
            >>> for kmer, tf in index.iter_kmers_by_frequency(min_tf=10, kmer_type="13mer"):
            ...     print(f"{kmer}: {tf}")
            
            >>> # Auto-detect mode and get top 1000
            >>> for kmer, tf in index.iter_kmers_by_frequency(max_kmers=1000):
            ...     print(f"{kmer}: {tf}")
        """
        if not self._loaded:
            raise RuntimeError("Index not loaded")
        
        # Auto-detect k-mer type if needed
        if kmer_type == "auto":
            # Check if 13-mer mode is available
            try:
                tf_array = self._wrapper.get_13mer_tf_array()
                kmer_type = "13mer"
            except:
                kmer_type = "23mer"
        
        if kmer_type == "13mer":
            # Use 13-mer specific method
            tf_array = self._wrapper.get_13mer_tf_array()
            
            # Create list of (index, tf) pairs for non-zero frequencies
            freq_list = []
            for index, tf in enumerate(tf_array):
                if tf >= min_tf:
                    freq_list.append((index, tf))
            
            # Sort by frequency (descending)
            freq_list.sort(key=lambda x: x[1], reverse=True)
            
            # Apply max_kmers limit if specified
            if max_kmers is not None:
                freq_list = freq_list[:max_kmers]
            
            # Yield k-mers with their frequencies
            for index, tf in freq_list:
                kmer = self._index_to_13mer(index)
                yield kmer, tf
                
        elif kmer_type == "23mer":
            # Use traditional method for 23-mers
            if not hasattr(self, 'n_kmers') or self.n_kmers == 0:
                raise RuntimeError("23-mer index not properly loaded")
            
            # Create list of (kmer, tf) pairs
            freq_list = []
            for kid in range(self.n_kmers):
                try:
                    kmer = self.get_kmer_by_kid(kid)
                    tf = self.get_tf_value(kmer)
                    if tf >= min_tf:
                        freq_list.append((kmer, tf))
                except:
                    continue
            
            # Sort by frequency (descending)
            freq_list.sort(key=lambda x: x[1], reverse=True)
            
            # Apply max_kmers limit if specified
            if max_kmers is not None:
                freq_list = freq_list[:max_kmers]
            
            # Yield k-mers with their frequencies
            for kmer, tf in freq_list:
                yield kmer, tf
        else:
            raise ValueError(f"Unsupported kmer_type: {kmer_type}. Use '13mer', '23mer', or 'auto'")

    def get_top_kmers(self, n: int = 100, min_tf: int = 1, kmer_type: str = "auto") -> List[Tuple[str, int]]:
        """
        Get top N most frequent k-mers
        
        Args:
            n: Number of top k-mers to return (default: 100)
            min_tf: Minimum term frequency threshold (default: 1)
            kmer_type: Type of k-mers: "13mer", "23mer", or "auto" (default: "auto")
            
        Returns:
            List of (kmer, term_frequency) tuples sorted by decreasing frequency
            
        Examples:
            >>> # Get top 50 most frequent 13-mers
            >>> top_kmers = index.get_top_kmers(50, kmer_type="13mer")
            >>> for kmer, tf in top_kmers:
            ...     print(f"{kmer}: {tf}")
        """
        return list(self.iter_kmers_by_frequency(min_tf=min_tf, max_kmers=n, kmer_type=kmer_type))

    def get_kmer_frequency_stats(self, kmer_type: str = "auto") -> Dict[str, Any]:
        """
        Get frequency statistics for k-mers
        
        Args:
            kmer_type: Type of k-mers: "13mer", "23mer", or "auto" (default: "auto")
            
        Returns:
            Dictionary with statistics: total_kmers, non_zero_kmers, max_tf, min_tf, avg_tf, etc.
            
        Examples:
            >>> stats = index.get_kmer_frequency_stats("13mer")
            >>> print(f"Total k-mers: {stats['total_kmers']}")
            >>> print(f"Non-zero k-mers: {stats['non_zero_kmers']}")
            >>> print(f"Max frequency: {stats['max_tf']}")
        """
        if not self._loaded:
            raise RuntimeError("Index not loaded")
        
        # Auto-detect k-mer type if needed
        if kmer_type == "auto":
            try:
                tf_array = self._wrapper.get_13mer_tf_array()
                kmer_type = "13mer"
            except:
                kmer_type = "23mer"
        
        if kmer_type == "13mer":
            tf_array = self._wrapper.get_13mer_tf_array()
            
            total_kmers = len(tf_array)
            non_zero_tf = [tf for tf in tf_array if tf > 0]
            non_zero_kmers = len(non_zero_tf)
            
            if non_zero_kmers > 0:
                max_tf = max(non_zero_tf)
                min_tf = min(non_zero_tf)
                avg_tf = sum(non_zero_tf) / non_zero_kmers
                total_tf = sum(tf_array)
            else:
                max_tf = min_tf = avg_tf = total_tf = 0
            
            return {
                'kmer_type': '13mer',
                'total_kmers': total_kmers,
                'non_zero_kmers': non_zero_kmers,
                'zero_kmers': total_kmers - non_zero_kmers,
                'max_tf': max_tf,
                'min_tf': min_tf,
                'avg_tf': avg_tf,
                'total_tf': total_tf,
                'coverage': non_zero_kmers / total_kmers if total_kmers > 0 else 0
            }
        
        elif kmer_type == "23mer":
            if not hasattr(self, 'n_kmers') or self.n_kmers == 0:
                raise RuntimeError("23-mer index not properly loaded")
            
            frequencies = []
            for kid in range(self.n_kmers):
                try:
                    kmer = self.get_kmer_by_kid(kid)
                    tf = self.get_tf_value(kmer)
                    frequencies.append(tf)
                except:
                    continue
            
            total_kmers = len(frequencies)
            non_zero_tf = [tf for tf in frequencies if tf > 0]
            non_zero_kmers = len(non_zero_tf)
            
            if non_zero_kmers > 0:
                max_tf = max(non_zero_tf)
                min_tf = min(non_zero_tf)
                avg_tf = sum(non_zero_tf) / non_zero_kmers
                total_tf = sum(frequencies)
            else:
                max_tf = min_tf = avg_tf = total_tf = 0
            
            return {
                'kmer_type': '23mer',
                'total_kmers': total_kmers,
                'non_zero_kmers': non_zero_kmers,
                'zero_kmers': total_kmers - non_zero_kmers,
                'max_tf': max_tf,
                'min_tf': min_tf,
                'avg_tf': avg_tf,
                'total_tf': total_tf,
                'coverage': non_zero_kmers / total_kmers if total_kmers > 0 else 0
            }
        else:
            raise ValueError(f"Unsupported kmer_type: {kmer_type}. Use '13mer', '23mer', or 'auto'")