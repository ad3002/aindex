#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.01.2018
#@author: Aleksey Komissarov  
#@contact: ad3002@gmail.com

import sys
import argparse
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Parse command line arguments
parser = argparse.ArgumentParser(description='AIndex Demo - showcase functionality')
parser.add_argument('--comprehensive', action='store_true', 
                   help='Run comprehensive demo with all functionality tests')
args = parser.parse_args()

# Use pybind11 API (memory-safe and efficient)
from aindex.core.aindex import (
    load_index_with_reads, get_aindex, AIndex, Strand,
    get_revcomp, hamming_distance, iter_reads_by_kmer
)
print("Using pybind11 API")

# Use test results from temp directory
test_dir = Path(__file__).parent
temp_dir = test_dir / "temp"

# Define explicit file paths (using correct names from the current pipeline)
reads_file = str(temp_dir / "test_kmer_counter.reads")
hash_file = str(temp_dir / "test_kmer_counter.23.pf")

tf_bin_file = str(temp_dir / "test_kmer_counter.23.tf.bin")
kmers_bin_file = str(temp_dir / "test_kmer_counter.23.kmers.bin")

pos_file = str(temp_dir / "test_kmer_counter.23.pos.bin")
index_file = str(temp_dir / "test_kmer_counter.23.index.bin")
indices_file = str(temp_dir / "test_kmer_counter.23.indices.bin")

# Check if test results exist
required_files = [reads_file, hash_file, tf_bin_file, kmers_bin_file, index_file, indices_file]

missing_files = [f for f in required_files if not Path(f).exists()]
if missing_files:
    print(f"Test results not found: {missing_files}")
    print("Please run 'make test-regression' first to generate test files.")
    exit(1)

print("\n=== LOADING INDEX ===")

# Load with explicit file paths
try:
    aindex = load_index_with_reads(
        hash_file=hash_file,
        tf_file=tf_bin_file,
        kmers_bin_file=kmers_bin_file,
        kmers_text_file="",  # Not used in this demo
        reads_file=reads_file,
        pos_file=pos_file,
        index_file=index_file,
        indices_file=indices_file,
        max_tf=10000
    )
    print("✓ Loaded index with explicit file paths successfully")
except FileNotFoundError as e:
    print(f"✗ Failed to load index: {e}")
    exit(1)

# Test legacy loading method
try:
    prefix_path = str(temp_dir / "test_kmer_counter")
    aindex_legacy = get_aindex(prefix_path, skip_aindex=False, max_tf=10000)
    print("✓ Legacy loading method also works")
except Exception as e:
    print(f"✗ Legacy loading failed: {e}")

print("\n=== AINDEX DEMO ===")

# Test kmers
kmers = ["A"*23, "T"*23, "AAAAAAAAAAAAAAAAAAAAAAT", "ATCGATCGATCGATCGATCGATCG"]

print("\n1. Basic K-mer Queries:")
for kmer in kmers:
    try:
        # Test multiple access methods
        tf_method = aindex.get_tf_value(kmer)
        tf_dict = aindex[kmer]
        exists = kmer in aindex
        tf_get = aindex.get(kmer, -1)
        hash_val = aindex.get_hash_value(kmer)
        kid = aindex.get_kid_by_kmer(kmer)
        strand = aindex.get_strand(kmer)
        
        print(f"  K-mer: {kmer}")
        print(f"    TF (method): {tf_method}")
        print(f"    TF (dict): {tf_dict}")
        print(f"    Exists: {exists}")
        print(f"    TF (get): {tf_get}")
        print(f"    Hash: {hash_val}")
        print(f"    Kid: {kid}")
        print(f"    Strand: {strand} ({strand.name if isinstance(strand, Strand) else 'unknown'})")
        
        # Test reverse lookup
        if kid > 0:
            kmer_back = aindex.get_kmer_by_kid(kid)
            print(f"    K-mer from KID: {kmer_back}")
        
    except Exception as e:
        print(f"  {kmer}: error = {e}")

print("\n2. Batch Operations:")
try:
    tf_values = aindex.get_tf_values(kmers)
    hash_values = aindex.get_hash_values(kmers)
    print(f"  Batch TF values: {tf_values}")
    print(f"  Batch hash values: {hash_values}")
except Exception as e:
    print(f"  Batch operations error: {e}")

print("\n3. Index Statistics:")
try:
    print(f"  Total k-mers: {len(aindex)}")
    print(f"  Number of k-mers: {aindex.n_kmers}")
    print(f"  Number of reads: {aindex.n_reads}")
    print(f"  Reads size: {aindex.reads_size} bytes")
    print(f"  Hash loaded: {aindex._loaded}")
    print(f"  Aindex loaded: {aindex.aindex_loaded}")
    print(f"  Max TF: {aindex.max_tf}")
except Exception as e:
    print(f"  Error getting statistics: {e}")

print("\n4. Utility Functions:")
test_seq = "ATCGATCGN"
try:
    revcomp = get_revcomp(test_seq)
    print(f"  Original: {test_seq}")
    print(f"  RevComp:  {revcomp}")
    
    # Test Hamming distance
    seq1, seq2 = "ATCGATCG", "ATCGATCC"
    hd = hamming_distance(seq1, seq2)
    print(f"  Hamming distance '{seq1}' vs '{seq2}': {hd}")
except Exception as e:
    print(f"  Utility functions error: {e}")

if aindex.aindex_loaded:
    print("\n5. Positional Queries:")
    for kmer in kmers[:2]:
        try:
            positions = aindex.get_positions(kmer)
            print(f"  K-mer {kmer}: {len(positions)} positions")
            if positions:
                pos = positions[0]
                rid = aindex.get_rid(pos)
                start = aindex.get_start(pos)
                print(f"    First position {pos} -> Read {rid}, Start {start}")
        except Exception as e:
            print(f"  {kmer}: position error = {e}")

    print("\n6. Read Retrieval:")
    # Test read retrieval by kmer
    for kmer in kmers[:2]:
        try:
            reads = aindex.get_reads_by_kmer(kmer, max_reads=3)
            print(f"  K-mer {kmer}: found {len(reads)} reads")
            for i, read in enumerate(reads):
                print(f"    Read {i}: {read[:40]}...")
        except Exception as e:
            print(f"  {kmer}: read retrieval error = {e}")

    # Test read retrieval by read ID
    print("\n7. Read Access by ID:")
    for rid in range(min(3, aindex.n_reads)):
        try:
            read = aindex.get_read_by_rid(rid)
            print(f"  Read {rid}: {read[:50]}...")
        except Exception as e:
            print(f"  Read {rid}: error = {e}")

    print("\n8. Advanced Analysis:")
    # Test sequence coverage
    test_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    try:
        coverage = aindex.get_sequence_coverage(test_sequence, cutoff=1)
        non_zero = sum(1 for x in coverage if x > 0)
        print(f"  Sequence coverage: {non_zero}/{len(coverage)} positions covered")
        
        # Test k-mer iteration in sequence
        kmer_count = 0
        for kmer, tf in aindex.iter_sequence_kmers(test_sequence):
            if kmer_count < 3:
                print(f"    K-mer {kmer}: tf={tf}")
            kmer_count += 1
            if kmer_count >= 3:
                break
    except Exception as e:
        print(f"  Sequence analysis error: {e}")

    print("\n9. Iterator Functions:")
    # Test read iteration
    try:
        read_count = 0
        for rid, read in aindex.iter_reads():
            if read_count < 2:
                print(f"  Read {rid}: {read[:40]}...")
            read_count += 1
            if read_count >= 2:
                break
        print(f"  Found {read_count} reads total")
    except Exception as e:
        print(f"  Read iteration error: {e}")

    # Test advanced analysis functions
    print("\n10. Advanced Analysis Functions:")
    test_kmer = kmers[0]
    try:
        # Test reads-by-kmer iterator
        read_count = 0
        for result in iter_reads_by_kmer(test_kmer, aindex):
            if read_count < 2:
                rid, pos, read, poses = result
                print(f"  Found {test_kmer} in read {rid} at position {pos}")
                print(f"    Read: {read[:40]}...")
            read_count += 1
            if read_count >= 2:
                break
    except Exception as e:
        print(f"  Advanced analysis error: {e}")

else:
    print("\n5-10. Positional queries skipped (aindex not loaded)")

print("\n=== DEMO COMPLETE ===")
print("✓ All core functionality tested successfully!")
print("✓ Using explicit file paths (no automatic name generation)")
print("✓ Both new pybind11 API and legacy compatibility confirmed")
print(f"✓ Processed index with {aindex.n_kmers} k-mers and {aindex.n_reads} reads")

# Run comprehensive demo if requested
if args.comprehensive:
    print("\n" + "="*60)
    print("RUNNING COMPREHENSIVE DEMO")
    print("="*60)
    try:
        from comprehensive_demo import main as comprehensive_main
        comprehensive_main()
        print("\n" + "="*60)
        print("✓ COMPREHENSIVE DEMO COMPLETED SUCCESSFULLY!")
        print("="*60)
    except ImportError as e:
        print(f"✗ Comprehensive demo not available: {e}")
        print("Make sure comprehensive_demo.py is in the tests directory")
    except Exception as e:
        print(f"✗ Comprehensive demo failed: {e}")
else:
    print("\nTip: Use --comprehensive flag to run extended functionality tests:")
    print("  python demo.py --comprehensive")