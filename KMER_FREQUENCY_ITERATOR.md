# Iterating Over k-mers Sorted by Frequency

The AIndex library now supports efficient iteration over k-mers sorted by frequency (starting from the most frequent).

## New Functions

### 1. `iter_kmers_by_frequency(min_tf=1, max_kmers=None, kmer_type="auto")`

An iterator that returns k-mers sorted in descending order of frequency.

**Parameters:**
- `min_tf` (int): Minimum k-mer frequency (default: 1)
- `max_kmers` (int, optional): Maximum number of k-mers to return (default: None - all)
- `kmer_type` (str): Type of k-mers - "13mer", "23mer", or "auto" (default: "auto")

**Returns:**
- Iterator[Tuple[str, int]]: Iterator of tuples (k-mer, frequency)

**Usage examples:**

```python
# Get all k-mers with frequency >= 100
for kmer, tf in index.iter_kmers_by_frequency(min_tf=100):
    print(f"{kmer}: {tf}")

# Get top 1000 most frequent 13-mers
for kmer, tf in index.iter_kmers_by_frequency(max_kmers=1000, kmer_type="13mer"):
    print(f"{kmer}: {tf}")

# Automatic k-mer type detection
for kmer, tf in index.iter_kmers_by_frequency(min_tf=50, max_kmers=500):
    print(f"{kmer}: {tf}")
```

### 2. `get_top_kmers(n=100, min_tf=1, kmer_type="auto")`

Returns a list of the top-N most frequent k-mers.

**Parameters:**
- `n` (int): Number of k-mers to return (default: 100)
- `min_tf` (int): Minimum k-mer frequency (default: 1)
- `kmer_type` (str): Type of k-mers - "13mer", "23mer", or "auto" (default: "auto")

**Returns:**
- List[Tuple[str, int]]: List of tuples (k-mer, frequency) sorted by descending frequency

**Usage examples:**

```python
# Get top 50 most frequent k-mers
top_kmers = index.get_top_kmers(n=50)
for kmer, tf in top_kmers:
    print(f"{kmer}: {tf}")

# Get top 200 13-mers with frequency >= 10
top_13mers = index.get_top_kmers(n=200, min_tf=10, kmer_type="13mer")
```

### 3. `get_kmer_frequency_stats(kmer_type="auto")`

Returns detailed statistics on k-mer frequencies.

**Parameters:**
- `kmer_type` (str): Type of k-mers - "13mer", "23mer", or "auto" (default: "auto")

**Returns:**
- Dict[str, Any]: Dictionary with statistics

**Result fields:**
- `kmer_type`: Type of k-mers ("13mer" or "23mer")
- `total_kmers`: Total number of k-mers
- `non_zero_kmers`: Number of k-mers with non-zero frequency
- `zero_kmers`: Number of k-mers with zero frequency
- `max_tf`: Maximum frequency
- `min_tf`: Minimum frequency (among non-zero)
- `avg_tf`: Average frequency (among non-zero)
- `total_tf`: Sum of all frequencies
- `coverage`: Fraction of k-mers with non-zero frequency (non_zero_kmers / total_kmers)

**Usage example:**

```python
stats = index.get_kmer_frequency_stats()
print(f"Type: {stats['kmer_type']}")
print(f"Total k-mers: {stats['total_kmers']:,}")
print(f"Non-zero: {stats['non_zero_kmers']:,}")
print(f"Coverage: {stats['coverage']:.4f}")
print(f"Max frequency: {stats['max_tf']:,}")
print(f"Average frequency: {stats['avg_tf']:.2f}")
```

## Supported k-mer Types

### 13-mers
For 13-mers, direct indexing of all possible combinations is used (4^13 = 67,108,864 k-mers). This provides very fast performance but requires loading a special 13-mer index.

### 23-mers
For 23-mers, a traditional hash index is used. It works slower than 13-mers but supports arbitrary sets of k-mers.

### Automatic Detection
When using `kmer_type="auto"`, the library automatically detects the available index type and uses it.

## Practical Usage Examples

### Frequency Distribution Analysis

```python
# Collect statistics on frequency distribution
freq_distribution = {}
for kmer, tf in index.iter_kmers_by_frequency(min_tf=1):
    if tf in freq_distribution:
        freq_distribution[tf] += 1
    else:
        freq_distribution[tf] = 1

# Show the most common frequencies
sorted_freqs = sorted(freq_distribution.items(), key=lambda x: x[1], reverse=True)
for tf, count in sorted_freqs[:10]:
    print(f"Frequency {tf}: {count} k-mers")
```

### Finding Rare k-mers

```python
# Find k-mers with frequency from 1 to 5
rare_kmers = []
for kmer, tf in index.iter_kmers_by_frequency(min_tf=1):
    if tf <= 5:
        rare_kmers.append((kmer, tf))
    if len(rare_kmers) >= 1000:  # Limit the number
        break

print(f"Found {len(rare_kmers)} rare k-mers")
```

### Exporting Results

```python
# Save top 10,000 k-mers to a file
with open('top_kmers.txt', 'w') as f:
    f.write("kmer\tfrequency\n")
    for kmer, tf in index.iter_kmers_by_frequency(max_kmers=10000):
        f.write(f"{kmer}\t{tf}\n")
```

### Analyzing High-Frequency k-mers

```python
# Find all k-mers with frequency above a threshold
threshold = 1000
high_freq_kmers = []

for kmer, tf in index.iter_kmers_by_frequency(min_tf=threshold):
    high_freq_kmers.append((kmer, tf))

print(f"K-mers with frequency >= {threshold}: {len(high_freq_kmers)}")

# Show top 20
for i, (kmer, tf) in enumerate(high_freq_kmers[:20]):
    print(f"{i+1:2d}. {kmer}: {tf:,}")
```

## Performance

- **13-mers**: Very fast iteration (~500,000 k-mers/sec) due to direct indexing
- **23-mers**: Moderate speed, depends on index size and number of unique k-mers

## Memory Requirements

- **13-mers**: Requires loading the full frequency array (~256 MB for uint32)
- **23-mers**: Memory usage depends on the number of unique k-mers in the index

## Compatibility

The new functions are fully compatible with the existing API and do not break old code.
