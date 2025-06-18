#!/usr/bin/env python3

"""
Quickstart: Iterating over k-mers sorted by frequency
"""

from aindex.core.aindex import AIndex

# 1. Load the index
index = AIndex()

# For 13-mers:
# index.load_13mer_index('path/to/hash.pf', 'path/to/frequencies.tf.bin')

# For 23-mers:
# index.load_hash('hash.pf', 'tf.bin', 'kmers.bin', 'kmers.txt')

# Automatic loading by prefix:
# index = AIndex.load_from_prefix('my_dataset')

# 2. Main operations

# Get top-100 most frequent k-mers
top_kmers = index.get_top_kmers(n=100)
print("Top-10 most frequent k-mers:")
for i, (kmer, tf) in enumerate(top_kmers[:10]):
    print(f"{i+1:2d}. {kmer}: {tf:,}")

# Iterate over k-mers with minimum frequency
print("\nK-mers with frequency >= 1000:")
count = 0
for kmer, tf in index.iter_kmers_by_frequency(min_tf=1000):
    print(f"{kmer}: {tf:,}")
    count += 1
    if count >= 5:  # Show only first 5
        break

# Get statistics
stats = index.get_kmer_frequency_stats()
print(f"\nIndex statistics:")
print(f"Type: {stats['kmer_type']}")
print(f"Total k-mers: {stats['total_kmers']:,}")
print(f"Non-zero k-mers: {stats['non_zero_kmers']:,}")
print(f"Max frequency: {stats['max_tf']:,}")
print(f"Average frequency: {stats['avg_tf']:.2f}")
print(f"Coverage: {stats['coverage']:.4f}")

# 3. Useful patterns

# Find rare k-mers (frequency 1-5)
print("\nRare k-mers (frequency 1-5):")
rare_count = 0
for kmer, tf in index.iter_kmers_by_frequency(min_tf=1):
    if tf <= 5:
        rare_count += 1
        if rare_count <= 10:  # Show first 10
            print(f"{kmer}: {tf}")
    if rare_count >= 100:  # Stop after 100
        break
print(f"Total {rare_count} rare k-mers found")

# Export to file
print("\nExporting top-1000 k-mers to file...")
with open('top_1000_kmers.txt', 'w') as f:
    f.write("kmer\tfrequency\n")
    for kmer, tf in index.iter_kmers_by_frequency(max_kmers=1000):
        f.write(f"{kmer}\t{tf}\n")
print("File 'top_1000_kmers.txt' created.")

print("\n" + "="*50)
print("Done! Now you can use new functions for k-mer analysis.")

