#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for analyzing k-mers in a file with reads.
Counts k-mer frequencies and their positions, creates a file with canonical k-mers.

@author: Aleksey Komissarov
@contact: ad3002@gmail.com
"""

import sys
import argparse
from collections import defaultdict
from typing import Dict, List, Tuple, Set


def reverse_complement(seq: str) -> str:
    """Get reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def get_canonical_kmer(kmer: str) -> str:
    """Get canonical form of k-mer (lexicographically smaller of forward and reverse)."""
    rev_comp = reverse_complement(kmer)
    return min(kmer, rev_comp)


def is_valid_kmer(kmer: str) -> bool:
    """Check if k-mer contains only valid nucleotides (A, T, G, C)."""
    return all(base in 'ATGC' for base in kmer)


def parse_fastq(filename: str) -> List[str]:
    """Parse FASTQ file and extract sequences."""
    sequences = []
    
    with open(filename, 'r') as f:
        line_count = 0
        for line in f:
            line = line.strip()
            if line_count % 4 == 1:  # Sequence line (2nd line in each block)
                sequences.append(line.upper())
            line_count += 1
    
    return sequences


def parse_multiple_fastq(filenames: List[str]) -> List[str]:
    """Parse multiple FASTQ files and combine sequences."""
    all_sequences = []
    
    for filename in filenames:
        print(f"Reading file: {filename}")
        sequences = parse_fastq(filename)
        all_sequences.extend(sequences)
        print(f"Read {len(sequences)} sequences from {filename}")
    
    return all_sequences


def analyze_kmers(sequences: List[str], k: int = 23) -> Tuple[Dict[str, int], Dict[str, List[Tuple[int, int, int]]]]:
    """
    Analyze k-mers in sequences.
    
    Args:
        sequences: List of DNA sequences
        k: K-mer length
    
    Returns:
        Tuple of:
        - Dictionary k-mer -> frequency
        - Dictionary k-mer -> list of positions (read_number, position, direction)
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
    Save results to file.
    
    Format: canonical_kmer frequency [read_number position direction] ...
    """
    with open(output_file, 'w') as f:
        # Sort k-mers by frequency (descending), then alphabetically
        sorted_kmers = sorted(kmer_counts.items(), 
                            key=lambda x: (-x[1], x[0]))
        
        for kmer, count in sorted_kmers:
            positions = kmer_positions[kmer]
            
            # Write k-mer and frequency
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
    parser = argparse.ArgumentParser(description='Analyze k-mers in file with reads')
    parser.add_argument('--input-files', nargs='+', 
                       default=['tests/data/exact100_reads_1.fasta', 'tests/data/exact100_reads_2.fasta'],
                       help='Input FASTQ files (default: exact100_reads_1.fasta and exact100_reads_2.fasta)')
    parser.add_argument('-k', '--kmer-size', type=int, default=23, 
                       help='K-mer size (default: 23)')
    parser.add_argument('-o', '--output', default='kmers_analysis.txt',
                       help='Output file (default: kmers_analysis.txt)')
    parser.add_argument('--min-count', type=int, default=1,
                       help='Minimum k-mer frequency for inclusion in results')
    
    args = parser.parse_args()
    
    print(f"Analyzing k-mers in files: {', '.join(args.input_files)}")
    print(f"K-mer size: {args.kmer_size}")
    print(f"Output file: {args.output}")
    
    # Parse FASTQ files
    print("Reading sequences...")
    sequences = parse_multiple_fastq(args.input_files)
    print(f"Total read {len(sequences)} sequences")
    
    # Analyze k-mers
    print("Analyzing k-mers...")
    kmer_counts, kmer_positions = analyze_kmers(sequences, args.kmer_size)
    
    # Filter by minimum frequency
    if args.min_count > 1:
        filtered_counts = {k: v for k, v in kmer_counts.items() if v >= args.min_count}
        filtered_positions = {k: v for k, v in kmer_positions.items() if k in filtered_counts}
        kmer_counts = filtered_counts
        kmer_positions = filtered_positions
        print(f"Applied minimum frequency filter: {args.min_count}")
    
    # Save results
    print("Saving results...")
    save_results(kmer_counts, kmer_positions, args.output)
    save_summary(kmer_counts, sequences, args.kmer_size, args.output)
    
    print(f"Analysis completed. Found {len(kmer_counts)} unique k-mers.")
    print(f"Results saved to: {args.output}")
    print(f"Summary saved to: {args.output}.summary")


if __name__ == "__main__":
    main()
