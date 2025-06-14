# AIndex Python API Documentation

## Overview

AIndex is a genomic data indexing library that provides efficient k-mer based indexing and querying of sequencing data. It uses perfect hashing for fast k-mer lookups and supports both hash-only and full positional indexing modes.

The library consists of:
- **C++ Core**: High-performance indexing algorithms implemented in C++
- **Python Bindings**: pybind11-based Python interface for easy integration
- **Legacy Support**: Backward-compatible API for existing codebases

## Architecture

### Core Components

1. **Perfect Hash Index**: Fast k-mer to frequency mapping using EMPHF library
2. **Positional Index**: Maps k-mers to their positions in reads
3. **Read Storage**: Memory-mapped read data with efficient access
4. **Interval Trees**: For mapping positions to read IDs and headers

### File Format

AIndex uses several file types:
- `.pf` - Perfect hash structure
- `.tf.bin` - Term frequency data (binary)
- `.kmers.bin` - K-mer binary data
- `.index.bin` - Position index
- `.indices.bin` - Position data
- `.pos.bin` - Position mappings
- `.reads` - Read sequences
- `.ridx` - Read index

## Python API Reference

### Main Classes

#### `AIndex`

The main class for genomic index operations.

```python
from aindex.core.aindex import AIndex

# Create instance
aindex = AIndex()
```

##### Constructor and Loading Methods

**`__init__()`**
Creates a new AIndex instance.

**`load_hash(hash_file: str, tf_file: str, kmers_bin_file: str, kmers_text_file: str)`**
Load hash index from explicit file paths.
- `hash_file`: Path to .pf file
- `tf_file`: Path to .tf.bin file
- `kmers_bin_file`: Path to .kmers.bin file
- `kmers_text_file`: Path to .kmers_text file (can be empty)

**`load_reads(reads_file: str)`**
Load reads file for positional queries.
- `reads_file`: Path to .reads file

**`load_aindex(pos_file: str, index_file: str, indices_file: str, max_tf: int)`**
Load positional index files.
- `pos_file`: Path to .pos.bin file
- `index_file`: Path to .index.bin file
- `indices_file`: Path to .indices.bin file
- `max_tf`: Maximum term frequency threshold

**`load_reads_index(index_file: str, header_file: str = None)`**
Load reads index and optional headers for interval mapping.
- `index_file`: Path to read index file
- `header_file`: Optional path to headers file

##### Core Query Methods

**`get_tf_value(kmer: str) -> int`**
Get term frequency for a single k-mer.
```python
tf = aindex.get_tf_value("ATCGATCGATCGATCGATCGATCG")
```

**`get_tf_values(kmers: List[str]) -> List[int]`**
Get term frequencies for multiple k-mers.
```python
tfs = aindex.get_tf_values(["ATCG...", "GCTA..."])
```

**`get_hash_value(kmer: str) -> int`**
Get hash value (k-mer ID) for a k-mer.

**`get_kid_by_kmer(kmer: str) -> int`**
Get k-mer ID by k-mer sequence.

**`get_kmer_by_kid(kid: int) -> str`**
Get k-mer sequence by k-mer ID.

**`get_strand(kmer: str) -> Strand`**
Get strand information for a k-mer.
Returns: `Strand.NOT_FOUND`, `Strand.FORWARD`, or `Strand.REVERSE`

##### Positional Query Methods

**`get_positions(kmer: str) -> List[int]`**
Get all positions where a k-mer occurs.
```python
positions = aindex.get_positions("ATCGATCGATCGATCGATCGATCG")
```

**`pos(kmer: str) -> List[int]`**
Alias for `get_positions()` (legacy compatibility).

**`get_rid(pos: int) -> int`**
Get read ID by position in reads file.

**`get_start(pos: int) -> int`**
Get start position of read by position.

**`get_read_by_rid(rid: int) -> str`**
Get read sequence by read ID.

**`get_read(start: int, end: int, revcomp: bool = False) -> str`**
Get read sequence by start and end positions.
- `revcomp`: If True, returns reverse complement

##### Advanced Query Methods

**`get_kmer_info(kid: int) -> Tuple[str, str, int]`**
Get comprehensive k-mer information by k-mer ID.
Returns: (kmer, reverse_complement_kmer, term_frequency)

**`get_rid2poses(kmer: str) -> Dict[int, List[int]]`**
Get mapping from read IDs to positions within reads for a k-mer.

**`get_reads_by_kmer(kmer: str, max_reads: int = 100) -> List[str]`**
Get reads containing a specific k-mer.

##### Iterator Methods

**`iter_reads()`**
Iterate over all reads.
```python
for rid, read in aindex.iter_reads():
    print(f"Read {rid}: {read}")
```

**`iter_reads_se()`**
Iterate over split-end reads (paired-end data).
```python
for rid, subread_idx, subread in aindex.iter_reads_se():
    print(f"Read {rid}.{subread_idx}: {subread}")
```

**`iter_sequence_kmers(sequence: str, k: int = 23)`**
Iterate over k-mers in a sequence with their frequencies.
```python
for kmer, tf in aindex.iter_sequence_kmers(sequence):
    print(f"{kmer}: {tf}")
```

##### Coverage Analysis Methods

**`get_sequence_coverage(seq: str, cutoff: int = 0, k: int = 23) -> List[int]`**
Get coverage array for a sequence based on k-mer frequencies.

**`print_sequence_coverage(seq: str, cutoff: int = 0)`**
Print sequence coverage and return coverage array.

##### Dictionary-like Interface

**`__len__() -> int`**
Get total number of k-mers in index.

**`__getitem__(kmer: str) -> int`**
Get term frequency for k-mer (dict-like access).
```python
tf = aindex["ATCGATCGATCGATCGATCGATCG"]
```

**`__contains__(kmer: str) -> bool`**
Check if k-mer exists in index.
```python
if "ATCGATCGATCGATCGATCGATCG" in aindex:
    print("K-mer found")
```

**`get(kmer: str, default: int = 0) -> int`**
Get term frequency with default value.

##### Properties

**`n_reads -> int`**
Number of reads in the index.

**`n_kmers -> int`**
Number of k-mers in the index.

**`aindex_loaded -> bool`**
Whether positional index is loaded.

**`reads_size -> int`**
Size of reads data.

### Utility Functions

#### `get_revcomp(sequence: str) -> str`
Get reverse complement of a DNA sequence.
```python
revcomp = get_revcomp("ATCG")  # Returns "CGAT"
```

#### `hamming_distance(s1: str, s2: str) -> int`
Calculate Hamming distance between two strings, ignoring 'N' positions.

### Loading Functions

#### `load_index(hash_file: str, tf_file: str, kmers_bin_file: str, kmers_text_file: str) -> AIndex`
Create and load an AIndex with explicit file paths.

#### `load_index_with_reads(hash_file: str, tf_file: str, ...) -> AIndex`
Create and load a fully functional AIndex with all components.

#### `get_aindex(prefix_path: str, skip_aindex: bool = False, max_tf: int = 1_000_000) -> AIndex`
Legacy function to load AIndex using filename prefix.
```python
# Loads files: prefix.23.pf, prefix.23.tf.bin, prefix.23.kmers.bin, etc.
aindex = get_aindex("data/my_index")
```

#### `load_aindex(settings: dict, ...) -> AIndex`
Flexible AIndex loading with configuration dictionary.

### Analysis Functions

#### `get_srandness(kmer: str, aindex: AIndex, k: int = 23) -> Tuple[int, int, int]`
Get strand distribution for a k-mer.
Returns: (plus_strand_count, minus_strand_count, total_count)

#### `iter_reads_by_kmer(kmer: str, aindex: AIndex, ...) -> Iterator`
Iterate over reads containing a specific k-mer.
```python
for rid, pos, read, poses in iter_reads_by_kmer(kmer, aindex):
    print(f"Found {kmer} in read {rid} at position {pos}")
```

#### `iter_reads_by_sequence(sequence: str, aindex: AIndex, ...) -> Iterator`
Iterate over reads containing a sequence with optional edit/Hamming distance.

#### `iter_reads_se_by_kmer(kmer: str, aindex: AIndex, ...) -> Iterator`
Iterate over split-end reads containing a k-mer.

#### `get_left_right_distances(left_kmer: str, right_kmer: str, aindex: AIndex, k: int = 23) -> Iterator`
Find distances between two k-mers in reads.

#### `get_layout_from_reads(kmer: str, aindex: AIndex, ...) -> Tuple`
Get alignment layout for reads containing a k-mer.

### Enumerations

#### `Strand`
Enumeration for strand information:
- `Strand.NOT_FOUND = 0`
- `Strand.FORWARD = 1`
- `Strand.REVERSE = 2`

## C++ Backend Reference

### `AindexWrapper` Class

The C++ backend class exposed through pybind11.

#### Core Methods

**`load(hash_file, tf_file, kmers_bin_file, kmers_text_file)`**
Load hash index files.

**`load_reads(reads_file)`**
Load reads file with memory mapping.

**`load_aindex(pos_file, index_file, indices_file, max_tf)`**
Load positional index files.

**`get_tf_value(kmer) -> uint64_t`**
Get term frequency for k-mer with strand checking.

**`get_kid_by_kmer(kmer) -> uint64_t`**
Get k-mer ID with reverse complement checking.

**`get_kmer_by_kid(kid) -> string`**
Convert k-mer ID back to sequence.

**`get_strand(kmer) -> uint64_t`**
Determine strand (0=not found, 1=forward, 2=reverse).

**`get_positions(kmer) -> vector<uint64_t>`**
Get all positions for k-mer up to max_tf limit.

#### Read Access Methods

**`get_read_by_rid(rid) -> string`**
Get read by read ID using start_positions mapping.

**`get_read(start, end, revcomp=false) -> string`**
Get read slice with optional reverse complement.

**`get_rid(pos) -> uint64_t`**
Map position to read ID using interval tree.

**`get_start(pos) -> uint64_t`**
Map position to read start position.

#### Data Structure Management

**`IntervalTree`**
Custom interval tree for position-to-read mapping:
- `addInterval(rid, start, end)`
- `query(start, end) -> vector<Interval>`

**`UsedReads`**
Set-based tracking of used read IDs with capacity limits.

#### Memory Management

- Memory-mapped file access for large data files
- Automatic cleanup in destructor
- Efficient interval tree queries

## Usage Examples

### Basic K-mer Frequency Lookup

```python
from aindex.core.aindex import get_aindex

# Load index
aindex = get_aindex("data/index_prefix")

# Query k-mer frequency
kmer = "ATCGATCGATCGATCGATCGATCG"
frequency = aindex[kmer]
print(f"K-mer {kmer} appears {frequency} times")

# Check if k-mer exists
if kmer in aindex:
    print("K-mer found in index")
```

### Positional Analysis

```python
# Get all positions of a k-mer
positions = aindex.get_positions(kmer)
print(f"K-mer found at {len(positions)} positions")

# Get reads containing the k-mer
for rid, pos, read, poses in iter_reads_by_kmer(kmer, aindex):
    print(f"Read {rid}: {read}")
    print(f"K-mer positions in read: {poses}")
```

### Coverage Analysis

```python
# Analyze sequence coverage
sequence = "ATCGATCGATCGATCGATCGATCGATCGATCG"
coverage = aindex.get_sequence_coverage(sequence, cutoff=5)

# Print coverage visualization
aindex.print_sequence_coverage(sequence, cutoff=5)
```

### Advanced Analysis

```python
# Find reads with specific sequence allowing mismatches
for result in iter_reads_by_sequence(target_seq, aindex, hd=2):
    rid, pos, read, poses, distance = result
    print(f"Found with {distance} mismatches in read {rid}")

# Analyze strand bias
plus, minus, total = get_srandness(kmer, aindex)
print(f"Strand distribution: +{plus}, -{minus}, total={total}")
```

## Performance Considerations

### Memory Usage
- Reads are memory-mapped for efficient access
- Perfect hash provides O(1) k-mer lookups
- Interval tree enables fast position-to-read mapping

### Scalability
- Handles genome-scale datasets (billions of k-mers)
- Configurable term frequency limits
- Efficient batch operations

### Best Practices
- Use `get_tf_values()` for batch k-mer queries
- Set appropriate `max_tf` limits for large datasets
- Use memory mapping for large read files
- Leverage iterator methods for memory-efficient processing

## File Dependencies

When loading an AIndex, ensure these files exist:
- `prefix.23.pf` - Perfect hash structure
- `prefix.23.tf.bin` - Term frequencies
- `prefix.23.kmers.bin` - K-mer data
- `prefix.23.index.bin` - Position index (for positional queries)
- `prefix.23.indices.bin` - Position data
- `prefix.23.pos.bin` - Position mappings
- `prefix.reads` - Read sequences
- `prefix.ridx` - Read index

Use explicit file paths with the new API for better control and error handling.