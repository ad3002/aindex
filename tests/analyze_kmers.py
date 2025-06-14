#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for analyzing kmers in a file with reads.
Counts kmer frequencies and their positions, creates a file with canonical kmers.
Processes files where each line contains one DNA sequence.

@author: Aleksey Komissarov
@contact: ad3002@gmail.com
"""

import sys
import os
import argparse
from collections import defaultdict
from typing import Dict, List, Tuple, Set


def reverse_complement(seq: str) -> str:
    """Get the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def get_canonical_kmer(kmer: str) -> str:
    """Get the canonical form of a kmer (lexicographically smallest of direct and reverse)."""
    rev_comp = reverse_complement(kmer)
    return min(kmer, rev_comp)


def is_valid_kmer(kmer: str) -> bool:
    """Check if the kmer contains only valid nucleotides (A, T, G, C)."""
    return all(base in 'ATGC' for base in kmer)


def parse_reads_file(filename: str) -> List[str]:
    """Parse a file with reads, where each line represents one sequence."""
    sequences = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line:  # Skip empty lines
                sequences.append(line.upper())
    
    return sequences


def analyze_kmers(sequences: List[str], k: int = 23) -> Tuple[Dict[str, int], Dict[str, List[Tuple[int, int, int]]]]:
    """
    Analyze kmers in sequences.
    
    Args:
        sequences: List of DNA sequences
        k: Kmer length
    
    Returns:
        Tuple of:
        - Dictionary kmer -> frequency
        - Dictionary kmer -> list of positions (read_number, position, direction)
    """
    kmer_counts = defaultdict(int)
    kmer_positions = defaultdict(list)
    
    for read_id, seq in enumerate(sequences):
        # Forward direction
        for pos in range(len(seq) - k + 1):
            kmer = seq[pos:pos + k]
            
            if not is_valid_kmer(kmer):
                continue
                
            canonical_kmer = get_canonical_kmer(kmer)
            kmer_counts[canonical_kmer] += 1
            
            # Determine direction: 0 - forward, 1 - reverse
            direction = 0 if kmer == canonical_kmer else 1
            kmer_positions[canonical_kmer].append((read_id, pos, direction))
    
    return dict(kmer_counts), dict(kmer_positions)


def save_results(kmer_counts: Dict[str, int], 
                kmer_positions: Dict[str, List[Tuple[int, int, int]]], 
                output_file: str):
    """
    Save results to a file.
    
    Format: canonical_kmer frequency [read_number position direction] ...
    """
    with open(output_file, 'w') as f:
        # Sort kmers by frequency (descending), then alphabetically
        sorted_kmers = sorted(kmer_counts.items(), 
                            key=lambda x: (-x[1], x[0]))
        
        for kmer, count in sorted_kmers:
            positions = kmer_positions[kmer]
            
            # Write kmer and frequency
            f.write(f"{kmer}\t{count}")
            
            # Write all positions
            for read_id, pos, direction in positions:
                f.write(f"\t{read_id},{pos},{direction}")
            
            f.write("\n")


def save_summary(kmer_counts: Dict[str, int], sequences: List[str], k: int, output_file: str):
    """Save summary statistics."""
    total_kmers = sum(kmer_counts.values())
    unique_kmers = len(kmer_counts)
    singleton_kmers = sum(1 for count in kmer_counts.values() if count == 1)
    max_count = max(kmer_counts.values()) if kmer_counts else 0
    
    with open(f"{output_file}.summary", 'w') as f:
        f.write("=== K-mer Analysis Summary ===\n")
        f.write(f"Input sequences: {len(sequences)}\n")
        f.write(f"K-mer size: {k}\n")
        f.write(f"Total k-mers: {total_kmers}\n")
        f.write(f"Unique k-mers: {unique_kmers}\n")
        f.write(f"Singleton k-mers: {singleton_kmers}\n")
        f.write(f"Max k-mer frequency: {max_count}\n")
        
        if unique_kmers > 0:
            avg_freq = total_kmers / unique_kmers
            f.write(f"Average frequency: {avg_freq:.2f}\n")
            
            # Theoretical maximum for given k
            theoretical_max = 4 ** k
            coverage = 100.0 * unique_kmers / theoretical_max
            f.write(f"K-mer space coverage: {coverage:.6f}%\n")


def main():
    parser = argparse.ArgumentParser(description='K-mer analysis in a file with reads')
    parser.add_argument('--input-file', 
                       default='./temp/reads.reads',
                       help='Input file with sequences (default: ./temp/reads.reads)')
    parser.add_argument('-k', '--kmer-size', type=int, default=23, 
                       help='K-mer size (default: 23)')
    parser.add_argument('-o', '--output', default='kmers_analysis.trues',
                       help='Output file (default: kmers_analysis.trues)')
    parser.add_argument('--min-count', type=int, default=1,
                       help='Minimum k-mer frequency to include in result')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: file {args.input_file} not found.")
        sys.exit(1)
    
    print(f"Analyzing kmers in file: {args.input_file}")
    print(f"K-mer size: {args.kmer_size}")
    print(f"Output file: {args.output}")
    
    # Parse reads file
    print("Reading sequences...")
    sequences = parse_reads_file(args.input_file)
    print(f"Total sequences read: {len(sequences)}")
    
    # Analyze kmers
    print("Analyzing kmers...")
    kmer_counts, kmer_positions = analyze_kmers(sequences, args.kmer_size)
    
    # Filter by minimum frequency
    if args.min_count > 1:
        filtered_counts = {k: v for k, v in kmer_counts.items() if v >= args.min_count}
        filtered_positions = {k: v for k, v in kmer_positions.items() if k in filtered_counts}
        kmer_counts = filtered_counts
        kmer_positions = filtered_positions
        print(f"Minimum frequency filter applied: {args.min_count}")
    
    # Save results
    print("Saving results...")
    save_results(kmer_counts, kmer_positions, args.output)
    save_summary(kmer_counts, sequences, args.kmer_size, args.output)
    
    print(f"Analysis complete. Found {len(kmer_counts)} unique kmers.")
    print(f"Results saved to: {args.output}")
    print(f"Summary saved to: {args.output}.summary")


if __name__ == "__main__":
    main()
