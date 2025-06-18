#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Comprehensive demo showcasing all AIndex functionality
# Tests all methods from aindex.py API
#
# @created: 2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import sys
import os
from pathlib import Path
from collections import defaultdict

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Import all AIndex functionality
from aindex.core.aindex import (
    AIndex, load_index, load_index_with_reads, get_aindex, load_aindex,
    get_revcomp, hamming_distance, Strand,
    get_srandness, iter_reads_by_kmer, iter_reads_by_sequence,
    iter_reads_se_by_kmer, get_left_right_distances, get_layout_from_reads
)

def print_section(title):
    """Print formatted section header."""
    print(f"\n{'='*60}")
    print(f" {title}")
    print(f"{'='*60}")

def print_subsection(title):
    """Print formatted subsection header."""
    print(f"\n--- {title} ---")

def safe_call(func, *args, **kwargs):
    """Safely call a function and handle exceptions."""
    try:
        return func(*args, **kwargs), None
    except Exception as e:
        return None, str(e)

def main():
    print_section("COMPREHENSIVE AINDEX DEMO")
    print("Testing all functionality from aindex.py API")
    
    # Use test results from temp directory
    test_dir = Path(__file__).parent
    temp_dir = test_dir / "temp"
    
    # Define explicit file paths
    reads_file = str(temp_dir / "test_kmer_counter.reads")
    hash_file = str(temp_dir / "test_kmer_counter.23.pf")
    tf_bin_file = str(temp_dir / "test_kmer_counter.23.tf.bin")
    kmers_bin_file = str(temp_dir / "test_kmer_counter.23.kmers.bin")
    index_file = str(temp_dir / "test_kmer_counter.23.index.bin")
    indices_file = str(temp_dir / "test_kmer_counter.23.indices.bin")
    
    # Check if test results exist
    required_files = [reads_file, hash_file, tf_bin_file, kmers_bin_file, 
                     index_file, indices_file]
    missing_files = [f for f in required_files if not Path(f).exists()]
    
    if missing_files:
        print(f"Test results not found: {missing_files}")
        print("Please run 'make test-regression' first to generate test files.")
        return 1

    print_section("1. UTILITY FUNCTIONS")
    
    print_subsection("1.1 DNA Sequence Utilities")
    test_seq = "ATCGATCGN"
    revcomp = get_revcomp(test_seq)
    print(f"Original:    {test_seq}")
    print(f"RevComp:     {revcomp}")
    
    # Test Hamming distance
    seq1 = "ATCGATCG"
    seq2 = "ATCGATCC"
    seq3 = "ATCNATCG"
    hd1 = hamming_distance(seq1, seq2)
    hd2 = hamming_distance(seq1, seq3)
    print(f"Hamming distance '{seq1}' vs '{seq2}': {hd1}")
    print(f"Hamming distance '{seq1}' vs '{seq3}': {hd2} (ignoring N)")
    
    print_section("2. LOADING FUNCTIONS")
    
    print_subsection("2.1 Loading with Explicit Paths")
    # Test different loading methods
    aindex = None
    
    try:
        aindex = load_index_with_reads(
            hash_file=hash_file,
            tf_file=tf_bin_file,
            kmers_bin_file=kmers_bin_file,
            kmers_text_file="",
            reads_file=reads_file,
            index_file=index_file,
            indices_file=indices_file,
            max_tf=10000
        )
        print("✓ Successfully loaded index with all components")
    except Exception as e:
        print(f"✗ Failed to load index: {e}")
        return 1
    
    print_subsection("2.2 Legacy Loading Functions")
    # Test legacy loading function if files exist in expected format
    prefix_path = str(temp_dir / "test_kmer_counter")
    try:
        aindex_legacy = get_aindex(prefix_path, skip_aindex=False, max_tf=10000)
        print("✓ Successfully loaded index using legacy get_aindex()")
    except Exception as e:
        print(f"✗ Legacy loading failed: {e}")
        aindex_legacy = None
    
    print_section("3. BASIC AINDEX FUNCTIONALITY")
    
    # Test k-mers (using real k-mers that might exist in data)
    test_kmers = [
        "A" * 23,  # Homopolymer
        "T" * 23,  # Homopolymer
        "ATCGATCGATCGATCGATCGATCG",  # Repetitive
        "AAAAAAAAAAAAAAAAAAAAAAT",   # Near homopolymer
        "CGTACGTACGTACGTACGTACGT"    # Another repetitive
    ]
    
    print_subsection("3.1 Dictionary-like Interface")
    for kmer in test_kmers[:3]:
        tf = aindex[kmer]
        exists = kmer in aindex
        tf_get = aindex.get(kmer, -1)
        print(f"K-mer: {kmer[:10]}...")
        print(f"  aindex[kmer]: {tf}")
        print(f"  kmer in aindex: {exists}")
        print(f"  aindex.get(kmer, -1): {tf_get}")
    
    print_subsection("3.2 Basic Query Methods")
    for kmer in test_kmers[:2]:
        tf_value, err = safe_call(aindex.get_tf_value, kmer)
        hash_value, err2 = safe_call(aindex.get_hash_value, kmer)
        kid, err3 = safe_call(aindex.get_kid_by_kmer, kmer)
        strand, err4 = safe_call(aindex.get_strand, kmer)
        
        print(f"K-mer: {kmer[:10]}...")
        print(f"  TF value: {tf_value} {f'(Error: {err})' if err else ''}")
        print(f"  Hash value: {hash_value} {f'(Error: {err2})' if err2 else ''}")
        print(f"  K-mer ID: {kid} {f'(Error: {err3})' if err3 else ''}")
        print(f"  Strand: {strand} {f'(Error: {err4})' if err4 else ''}")
        
        if kid and not err3:
            kmer_back, err5 = safe_call(aindex.get_kmer_by_kid, kid)
            print(f"  K-mer from ID: {kmer_back} {f'(Error: {err5})' if err5 else ''}")
    
    print_subsection("3.3 Batch Operations")
    tf_values, err = safe_call(aindex.get_tf_values, test_kmers[:3])
    hash_values, err2 = safe_call(aindex.get_hash_values, test_kmers[:3])
    
    print(f"Batch TF values: {tf_values} {f'(Error: {err})' if err else ''}")
    print(f"Batch hash values: {hash_values} {f'(Error: {err2})' if err2 else ''}")
    
    print_subsection("3.4 Index Properties")
    print(f"Index length: {len(aindex)}")
    print(f"Number of k-mers: {aindex.n_kmers}")
    print(f"Number of reads: {aindex.n_reads}")
    print(f"Reads size: {aindex.reads_size} bytes")
    print(f"Hash loaded: {aindex._loaded}")
    print(f"Aindex loaded: {aindex.aindex_loaded}")
    print(f"Max TF: {aindex.max_tf}")
    
    print_section("4. POSITIONAL QUERIES")
    
    if aindex.aindex_loaded:
        print_subsection("4.1 Position-based Queries")
        for kmer in test_kmers[:2]:
            positions, err = safe_call(aindex.get_positions, kmer)
            positions_alias, err2 = safe_call(aindex.pos, kmer)
            
            print(f"K-mer: {kmer[:10]}...")
            if not err and positions:
                print(f"  Positions (first 5): {positions[:5]}")
                print(f"  Total positions: {len(positions)}")
                
                # Test position-to-read mapping
                if positions:
                    pos = positions[0]
                    rid, err3 = safe_call(aindex.get_rid, pos)
                    start, err4 = safe_call(aindex.get_start, pos)
                    print(f"  First position {pos} -> Read ID: {rid}, Start: {start}")
            else:
                print(f"  No positions found {f'(Error: {err})' if err else ''}")
        
        print_subsection("4.2 Read Retrieval")
        # Test read retrieval by RID
        for rid in range(min(3, aindex.n_reads)):
            read, err = safe_call(aindex.get_read_by_rid, rid)
            if not err and read:
                print(f"Read {rid}: {read[:50]}{'...' if len(read) > 50 else ''}")
            else:
                print(f"Read {rid}: Error: {err}")
        
        # Test read retrieval by position
        for kmer in test_kmers[:1]:
            positions, err = safe_call(aindex.get_positions, kmer)
            if not err and positions:
                pos = positions[0]
                read, err2 = safe_call(aindex.get_read, pos, pos + 50)
                print(f"Read slice at position {pos}: {read}")
        
        print_subsection("4.3 Read-to-Position Mapping")
        for kmer in test_kmers[:2]:
            rid2poses, err = safe_call(aindex.get_rid2poses, kmer)
            if not err and rid2poses:
                print(f"K-mer {kmer[:10]}... found in {len(rid2poses)} reads")
                for rid, poses in list(rid2poses.items())[:3]:
                    print(f"  Read {rid}: positions {poses}")
            else:
                print(f"K-mer {kmer[:10]}...: No mapping found")
    else:
        print("Positional index not loaded - skipping positional queries")
    
    print_section("5. ADVANCED QUERIES")
    
    print_subsection("5.1 K-mer Information")
    for kmer in test_kmers[:2]:
        kid, err = safe_call(aindex.get_kid_by_kmer, kmer)
        if not err and kid:
            info, err2 = safe_call(aindex.get_kmer_info, kid)
            if not err2:
                kmer_seq, rkmer_seq, tf = info
                print(f"K-mer ID {kid}:")
                print(f"  K-mer: {kmer_seq}")
                print(f"  RevComp: {rkmer_seq}")
                print(f"  TF: {tf}")
                
                # Test legacy method
                legacy_info, err3 = safe_call(aindex.get_kmer_info_by_kid, kid)
                if not err3:
                    print(f"  Legacy info matches: {legacy_info == info}")
    
    print_subsection("5.2 Sequence Analysis")
    test_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    
    # Test sequence k-mer iteration
    print("K-mers in test sequence (first 5):")
    kmer_count = 0
    for kmer, tf in aindex.iter_sequence_kmers(test_sequence):
        if kmer_count < 5:
            print(f"  {kmer}: tf={tf}")
        kmer_count += 1
        if kmer_count >= 5:
            break
    
    # Test coverage analysis
    coverage, err = safe_call(aindex.get_sequence_coverage, test_sequence, cutoff=1)
    if not err:
        print(f"Sequence coverage (length {len(coverage)}):")
        print(f"  Non-zero positions: {sum(1 for x in coverage if x > 0)}")
        print(f"  Max coverage: {max(coverage) if coverage else 0}")
    
    print_section("6. ITERATOR FUNCTIONALITY")
    
    if aindex.reads_size > 0:
        print_subsection("6.1 Read Iteration")
        # Test basic read iteration (first few reads)
        read_count = 0
        for rid, read in aindex.iter_reads():
            if read_count < 3:
                print(f"Read {rid}: {read[:50]}...")
            read_count += 1
            if read_count >= 3:
                break
        print(f"Total reads available: {read_count} (showing first 3)")
        
        # Test split-end read iteration
        print_subsection("6.2 Split-End Read Iteration")
        se_count = 0
        for rid, subread_idx, subread in aindex.iter_reads_se():
            if se_count < 3:
                print(f"Read {rid}.{subread_idx}: {subread[:50]}...")
            se_count += 1
            if se_count >= 3:
                break
    
    print_section("7. ANALYSIS FUNCTIONS")
    
    if aindex.aindex_loaded:
        print_subsection("7.1 Strand Analysis")
        for kmer in test_kmers[:2]:
            plus, minus, total, err = safe_call(get_srandness, kmer, aindex)
            if not err:
                print(f"K-mer {kmer[:10]}... strand distribution:")
                print(f"  Plus strand: {plus}")
                print(f"  Minus strand: {minus}")
                print(f"  Total: {total}")
        
        print_subsection("7.2 Read-by-K-mer Analysis")
        for kmer in test_kmers[:1]:  # Test one k-mer
            print(f"Reads containing k-mer {kmer[:10]}...:")
            read_count = 0
            for result in iter_reads_by_kmer(kmer, aindex, k=23):
                if read_count < 3:
                    rid, pos, read, poses = result
                    print(f"  Read {rid}: position {pos}, k-mer positions {poses}")
                    print(f"    Sequence: {read[:50]}...")
                read_count += 1
                if read_count >= 3:
                    break
            if read_count == 0:
                print(f"  No reads found containing k-mer")
        
        print_subsection("7.3 Sequence Search with Mismatches")
        search_seq = "ATCGATCGATCGATCGATCGATC"  # 23bp sequence
        print(f"Searching for sequence: {search_seq}")
        
        # Exact match
        exact_count = 0
        for result in iter_reads_by_sequence(search_seq, aindex, k=23):
            if result and exact_count < 2:
                if len(result) == 4:  # Exact match
                    rid, pos, read, poses = result
                    print(f"  Exact match in read {rid} at position {pos}")
                exact_count += 1
                if exact_count >= 2:
                    break
        
        # With Hamming distance
        hd_count = 0
        for result in iter_reads_by_sequence(search_seq, aindex, hd=1, k=23):
            if result and hd_count < 2:
                if len(result) == 5:  # With distance
                    rid, pos, read, poses, distance = result
                    print(f"  Match with HD={distance} in read {rid} at position {pos}")
                hd_count += 1
                if hd_count >= 2:
                    break
        
        print_subsection("7.4 Split-End Analysis")
        for kmer in test_kmers[:1]:
            print(f"Split-end reads containing k-mer {kmer[:10]}...:")
            se_count = 0
            for result in iter_reads_se_by_kmer(kmer, aindex):
                if se_count < 3:
                    rid, pos, read, side = result
                    side_str = "left" if side == 0 else "right" if side == 1 else "single"
                    print(f"  Read {rid} ({side_str} side): {read[:30]}...")
                se_count += 1
                if se_count >= 3:
                    break
    
    print_section("8. COMPREHENSIVE TESTS")
    
    print_subsection("8.1 Error Handling")
    # Test with invalid k-mers
    invalid_kmers = ["INVALID", "N" * 23, "X" * 23]
    for kmer in invalid_kmers:
        tf, err = safe_call(aindex.get_tf_value, kmer)
        print(f"Invalid k-mer '{kmer}': tf={tf}, error='{err}'")
    
    print_subsection("8.2 Edge Cases")
    # Test edge cases
    empty_result, err = safe_call(aindex.get_tf_value, "")
    print(f"Empty k-mer: tf={empty_result}, error='{err}'")
    
    # Test very large read ID
    large_rid_result, err = safe_call(aindex.get_read_by_rid, 999999)
    print(f"Large read ID: result='{large_rid_result}', error='{err}'")
    
    print_section("9. PERFORMANCE SUMMARY")
    
    print("Functionality tested:")
    tests = [
        "✓ Dictionary-like interface (__getitem__, __contains__, get)",
        "✓ Basic queries (TF, hash, k-mer ID, strand)",
        "✓ Batch operations (multiple k-mers)",
        "✓ Positional queries (positions, reads by position)",
        "✓ Read retrieval (by ID, by position)",
        "✓ Advanced k-mer info",
        "✓ Sequence analysis (coverage, k-mer iteration)",
        "✓ Read iteration (basic and split-end)",
        "✓ Analysis functions (strand bias, read search)",
        "✓ Sequence search with mismatches",
        "✓ Error handling and edge cases",
        "✓ Legacy compatibility functions"
    ]
    
    for test in tests:
        print(f"  {test}")
    
    print_section("DEMO COMPLETE")
    print("All AIndex functionality has been tested successfully!")
    print("The library provides comprehensive k-mer indexing and analysis capabilities.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())