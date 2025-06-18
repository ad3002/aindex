#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Example usage of k-mer frequency iteration functions
"""

from aindex.core.aindex import AIndex

def test_kmer_frequency_iteration():
    """
    Example usage of new functions for iterating over k-mers
    """
    
    # Load index (replace with real file paths)
    index = AIndex()
    
    # Example for 13-mers
    try:
        # Load 13-mer index
        # index.load_13mer_index('path/to/13mers.pf', 'path/to/13mers.tf.bin')
        
        print("=== Examples of frequency iterator usage ===")
        
        # 1. Get top-100 most frequent k-mers
        print("\n1. Top-100 most frequent k-mers:")
        top_kmers = index.get_top_kmers(n=100, kmer_type="13mer")
        for i, (kmer, tf) in enumerate(top_kmers[:10]):  # Show only first 10
            print(f"{i+1}. {kmer}: {tf}")
        
        # 2. Iterate over k-mers with minimum frequency
        print("\n2. K-mers with frequency >= 100:")
        count = 0
        for kmer, tf in index.iter_kmers_by_frequency(min_tf=100, kmer_type="13mer"):
            print(f"{kmer}: {tf}")
            count += 1
            if count >= 20:  # Limit output
                print("... and others")
                break
        
        # 3. Get frequency statistics
        print("\n3. K-mer frequency statistics:")
        stats = index.get_kmer_frequency_stats(kmer_type="13mer")
        print(f"K-mer type: {stats['kmer_type']}")
        print(f"Total k-mers: {stats['total_kmers']:,}")
        print(f"Non-zero k-mers: {stats['non_zero_kmers']:,}")
        print(f"Zero k-mers: {stats['zero_kmers']:,}")
        print(f"Max frequency: {stats['max_tf']:,}")
        print(f"Min frequency: {stats['min_tf']:,}")
        print(f"Average frequency: {stats['avg_tf']:.2f}")
        print(f"Total frequency: {stats['total_tf']:,}")
        print(f"Coverage: {stats['coverage']:.4f}")
        
        # 4. Automatic k-mer type detection
        print("\n4. Automatic k-mer type detection:")
        auto_top = index.get_top_kmers(n=10, kmer_type="auto")
        for kmer, tf in auto_top:
            print(f"{kmer}: {tf}")
            
    except Exception as e:
        print(f"Error during testing: {e}")
        print("Make sure the index is loaded correctly")

def example_usage_patterns():
    """
    Additional usage examples
    """
    
    print("\n=== Additional usage examples ===")
    
    # Example code for different scenarios
    examples = [
        {
            "title": "Finding rare k-mers",
            "code": """
# Find k-mers with frequency from 1 to 5
rare_kmers = []
for kmer, tf in index.iter_kmers_by_frequency(min_tf=1, kmer_type="13mer"):
    if tf <= 5:
        rare_kmers.append((kmer, tf))
    if len(rare_kmers) >= 100:
        break
print(f"Found {len(rare_kmers)} rare k-mers")
"""
        },
        {
            "title": "Frequency distribution analysis",
            "code": """
# Collect frequency distribution statistics
freq_distribution = {}
for kmer, tf in index.iter_kmers_by_frequency(min_tf=1, kmer_type="13mer"):
    if tf in freq_distribution:
        freq_distribution[tf] += 1
    else:
        freq_distribution[tf] = 1

# Show top-10 most common frequency values
sorted_freqs = sorted(freq_distribution.items(), key=lambda x: x[1], reverse=True)
for tf, count in sorted_freqs[:10]:
    print(f"Frequency {tf}: {count} k-mers")
"""
        },
        {
            "title": "Threshold comparison",
            "code": """
# Count k-mers above a threshold
threshold = 100
high_freq_count = 0
for kmer, tf in index.iter_kmers_by_frequency(min_tf=threshold, kmer_type="13mer"):
    high_freq_count += 1

print(f"K-mers with frequency >= {threshold}: {high_freq_count}")
"""
        },
        {
            "title": "Export to file",
            "code": """
# Save top-1000 k-mers to a file
with open('top_kmers.txt', 'w') as f:
    f.write("kmer\\tfrequency\\n")
    for kmer, tf in index.iter_kmers_by_frequency(max_kmers=1000, kmer_type="13mer"):
        f.write(f"{kmer}\\t{tf}\\n")
"""
        }
    ]
    
    for example in examples:
        print(f"\n{example['title']}:")
        print(example['code'])

if __name__ == "__main__":
    print("Demonstration of k-mer frequency iteration functions")
    print("=" * 60)
    
    # Main examples
    test_kmer_frequency_iteration()
    
    # Additional examples
    example_usage_patterns()
    
    print("\n" + "=" * 60)
    print("Documentation:")
    print("""
Main functions:

1. iter_kmers_by_frequency(min_tf=1, max_kmers=None, kmer_type="auto")
   - Iterator over k-mers sorted by descending frequency
   - min_tf: minimum frequency (default 1)
   - max_kmers: maximum number of k-mers (default None - all)
   - kmer_type: k-mer type ("13mer", "23mer", "auto")

2. get_top_kmers(n=100, min_tf=1, kmer_type="auto")
   - Returns a list of top-N most frequent k-mers
   - n: number of k-mers to return
   - min_tf: minimum frequency
   - kmer_type: k-mer type

3. get_kmer_frequency_stats(kmer_type="auto")
   - Returns k-mer frequency statistics
   - Includes: total count, non-zero, max/min/average frequencies, coverage

All functions support automatic k-mer type detection (13-mer or 23-mer).
""")
