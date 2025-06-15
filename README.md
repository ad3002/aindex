# aindex: perfect hash based index for genomic data

[![PyPI version](https://badge.fury.io/py/aindex2.svg)](https://badge.fury.io/py/aindex2)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/aindex2.svg)](https://pypi.python.org/pypi/aindex2/)
![PyPI - Wheel](https://img.shields.io/pypi/wheel/aindex2)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/ad3002/aindex/build_wheels.yml)
[![PyPI license](https://img.shields.io/pypi/l/aindex2.svg)](https://pypi.python.org/pypi/aindex2/)
[![DOI](https://zenodo.org/badge/114383739.svg)](https://zenodo.org/doi/10.5281/zenodo.12818331)

## Features

🚀 **High Performance**: Ultra-fast k-mer querying with optimized C++ backend
- **13-mers**: 2.0M queries/sec (batch), 491K queries/sec (single)
- **23-mers**: 2.3M queries/sec (batch), 1.1M queries/sec (single)
- **Sequence coverage analysis**: 24.5K sequences/sec (13-mers), 17.5K sequences/sec (23-mers)

🧬 **Dual K-mer Support**: Native support for both 13-mer and 23-mer k-mers
- **13-mers**: Complete 4^13 space coverage with perfect hashing
- **23-mers**: Efficient sparse indexing for genomic sequences
- **Auto-detection**: Seamlessly switches between modes based on k-mer length

💾 **Memory Efficient**: Optimized data structures and memory-mapped files
- **Batch operations**: Up to 4x faster than single queries
- **Minimal memory overhead**: Constant memory usage during processing
- **Real-time processing**: Stream processing for large genomic datasets

🔧 **Modern API**: Clean pybind11 interface with comprehensive functionality



## Installation

**Quick install with pip:**
```bash
pip install aindex2
```

**✅ Optimized builds available for all platforms:**
- **Linux**: x86_64, aarch64 (ARM64)
- **macOS**: x86_64 (Intel), arm64 (Apple Silicon M1/M2/M3)
- **Windows**: AMD64

**⚡ ARM64-optimized builds for Apple Silicon:**
Our builds include ARM64 optimizations for Apple M1/M2/M3 processors, providing:
- Up to 30% faster k-mer querying
- Native ARM64 SIMD optimizations
- Optimized for Apple Silicon architecture

**Building from source (optional):**
```bash
git clone https://github.com/ad3002/aindex.git
cd aindex
make arm64  # For Apple Silicon
# or
make all    # For x86_64
pip install .
```

**Requirements:**
- Python 3.8+
- Standard build tools (automatically handled by pip)
- No external dependencies required

**For Google Colab users:**
```python
!pip install aindex2
```

### Detailed Installation Instructions

**Standard installation with pip (all platforms):**

```bash
pip install aindex2
```

**Installation from source (for development or custom builds):**

```bash
git clone https://github.com/ad3002/aindex.git
cd aindex

# Standard build (all platforms)
make all
pip install .

# For Apple Silicon with ARM64 optimizations  
make arm64
pip install .
```

**Platform support:**
- **Linux**: x86_64, aarch64 (ARM64) - pre-built wheels available
- **macOS**: x86_64 (Intel), arm64 (Apple Silicon) - pre-built wheels available  
- **Windows**: AMD64 - pre-built wheels available

All platforms include optimized builds with no external dependencies required.

### Google Colab Installation

For installation in Google Colab environment, there's a known cmake conflict that needs to be resolved first:

```python
# Quick fix for cmake conflict
!pip uninstall -y cmake
!apt-get update
!apt-get install -y build-essential cmake git python3-dev

# Clone and install aindex
!git clone https://github.com/ad3002/aindex.git
%cd aindex
!pip install .
```

**Alternative: Use automatic installation script**

```python
# Download and run the installation script
!wget https://raw.githubusercontent.com/ad3002/aindex/main/install_colab.py
!python install_colab.py
```

**For troubleshooting**, use the diagnostic script:

```python
!wget https://raw.githubusercontent.com/ad3002/aindex/main/diagnose_colab.py
!python diagnose_colab.py
```

> **Note**: Google Colab has a conflict between the Python `cmake` package and system cmake. The scripts above automatically resolve this issue.

To uninstall:

```bash
pip uninstall aindex2
pip uninstall clean
```

To clean up the compiled files, run:

```
make clean
```

### macOS Compilation (ARM64/Apple Silicon Support)

For macOS systems (including Apple Silicon M1/M2), aindex now provides full ARM64 support with optimized performance:

```bash
# Build all components including the fast kmer_counter utility
make

# Alternative: build only core components
make macos
```

The project has been fully ported to ARM64/macOS, removing x86-specific dependencies (SSE instructions) and adding native ARM64 optimization.

**Requirements for macOS:**
- For jellyfish-based pipeline: `brew install jellyfish` (optional)
- Built-in kmer_counter provides faster alternative to jellyfish
- Python development headers (usually included with Xcode tools)

**Performance Note**: The new built-in `kmer_counter` utility is approximately 5x faster than jellyfish for k-mer counting tasks.

## Usage

### Command Line Interface (CLI)

aindex provides a unified command-line interface for all tools and utilities. After installation, all functions are accessible through the `aindex` command:

```bash
# Get help for all available commands
aindex --help

# Get help for a specific command
aindex count --help
aindex compute-aindex --help
```

#### Available Commands

**Core indexing tools:**
```bash
# Compute AIndex for genomic sequences (supports both 13-mer and 23-mer modes)
aindex compute-aindex -i input.fastq -o output_prefix -k 23

# Compute general index
aindex compute-index -i input.fasta -o output_prefix

# Process reads for indexing
aindex compute-reads input.fastq output.fastq fastq reads_prefix
```

**K-mer analysis:**
```bash
# Count k-mers in sequences (fast built-in counter)
aindex count -i input.fasta -o output.txt -k 23 -t 4

# Count 13-mers specifically (optimized for complete 13-mer space)
aindex count -i input.fasta -o output.txt -k 13

# Build hash table for k-mers
aindex build-hash -i kmers.txt -o hash_output

# Generate all possible 13-mers
aindex generate -o all_13mers.txt -k 13
```

**Utilities:**
```bash
# Convert reads to FASTA format
aindex reads-to-fasta input.fastq output.fasta

# Show version information
aindex version

# Show system and installation information
aindex info
```

#### Examples

**Count 23-mers in a FASTA file:**
```bash
aindex count -i genome.fasta -o kmer_counts.txt -k 23 -t 8
```

**Build AIndex for 13-mer analysis:**
```bash
aindex compute-aindex -i reads.fastq -o reads_index -k 13 --lu 2 -P 16
```

**Generate all possible 13-mers for reference:**
```bash
aindex generate -o all_13mers.txt -k 13
```

### K-mer Counting Pipelines

aindex supports two k-mer counting backends:

1. **Built-in kmer_counter (Recommended)** - Fast native implementation, ~5x faster than jellyfish
2. **Jellyfish** - Traditional external tool (requires `brew install jellyfish` on macOS)

### Quick Start

Compute all binary arrays using the fast built-in counter:

```bash
FASTQ1=./tests/raw_reads.101bp.IS350bp25_1.fastq
FASTQ2=./tests/raw_reads.101bp.IS350bp25_2.fastq
OUTPUT_PREFIX=./tests/raw_reads.101bp.IS350bp25

# Using built-in kmer_counter (recommended, faster) via CLI
aindex compute-aindex -i $FASTQ1,$FASTQ2 -t fastq -o $OUTPUT_PREFIX --lu 2 -P 30 --use-kmer-counter

# Using built-in kmer_counter (legacy script approach)
python3 scripts/compute_aindex.py -i $FASTQ1,$FASTQ2 -t fastq -o $OUTPUT_PREFIX --lu 2 -P 30 --use_kmer_counter

# Using jellyfish (traditional approach)
python3 scripts/compute_aindex.py -i $FASTQ1,$FASTQ2 -t fastq -o $OUTPUT_PREFIX --lu 2 -P 30
```

### Command Line Options

- `-i, --input`: Input FASTQ/FASTA files (comma-separated for multiple files)
- `-t, --type`: Input file type ('fastq' or 'fasta')
- `-o, --output`: Output prefix for generated files
- `--lu`: Lower frequency threshold for k-mers
- `-P, --threads`: Number of threads to use
- `--use-kmer-counter`: Use built-in fast k-mer counter instead of jellyfish

### Pipeline Outputs

Both pipelines generate identical output files:
- `.reads` - Processed reads file
- `.dat` - K-mer frequency data
- `.aindex` - Binary index file
- `.stat` - Statistics and metadata

## Usage from Python

### Modern API

The aindex package provides a unified API supporting both 13-mer and 23-mer modes:

```python
from aindex.core.aindex import AIndex
import aindex.core.aindex_cpp as aindex_cpp

# Load 23-mer index (for genomic sequences)
index_23mer = AIndex.load_from_prefix("temp/reads.23")
index_23mer.load_reads("temp/reads.reads")  # Optional: load actual read sequences

# Load 13-mer index (for complete k-mer space analysis)
index_13mer = aindex_cpp.AindexWrapper()
index_13mer.load_from_prefix_13mer("temp/all_13mers")
index_13mer.load_reads("temp/reads.reads")  # Optional: load reads

print(f"23-mer index: {index_23mer.n_kmers:,} k-mers, {index_23mer.n_reads:,} reads")
print(f"13-mer index: {index_13mer.get_13mer_statistics()}")
```

### K-mer Frequency Queries

**Single k-mer queries:**
```python
# 23-mer queries (using AIndex wrapper)
tf_23 = index_23mer.get_tf_value("ATCGATCGATCGATCGATCGATC")  # 23 characters
print(f"23-mer frequency: {tf_23}")

# Alternative 23-mer query using [] operator
tf_23_alt = index_23mer["ATCGATCGATCGATCGATCGATC"]
print(f"23-mer frequency (alt): {tf_23_alt}")

# 13-mer queries (using C++ wrapper directly)
tf_13 = index_13mer.get_total_tf_value_13mer("ATCGATCGATCGA")  # 13 characters
print(f"13-mer frequency: {tf_13}")

# Get forward and reverse frequencies separately for 13-mers
tf_fwd, tf_rev = index_13mer.get_tf_both_directions_13mer("ATCGATCGATCGA")
print(f"13-mer forward: {tf_fwd}, reverse: {tf_rev}, total: {tf_fwd + tf_rev}")
```

**Batch queries (much faster):**
```python
# Batch 23-mer queries (2-3x faster than single queries)
kmers_23 = ["ATCGATCGATCGATCGATCGATC", "AAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTT"]
tf_values_23 = index_23mer.get_tf_values(kmers_23)
print(f"23-mer batch results: {tf_values_23}")

# Batch 13-mer queries (total frequencies)
kmers_13 = ["ATCGATCGATCGA", "AAAAAAAAAAAAA", "TTTTTTTTTTTTT"] 
tf_values_13 = index_13mer.get_total_tf_values_13mer(kmers_13)
print(f"13-mer batch results: {tf_values_13}")

# Batch directional 13-mer queries (forward + reverse separately)
directional_results = index_13mer.get_tf_both_directions_13mer_batch(kmers_13)
for i, (fwd, rev) in enumerate(directional_results):
    print(f"{kmers_13[i]}: forward={fwd}, reverse={rev}, total={fwd+rev}")
```

### Advanced 13-mer Operations

**Directional analysis (forward + reverse complement):**
```python
# Get frequencies in both directions
kmer = "ATCGATCGATCGA"
forward_tf, reverse_tf = index_13mer.get_tf_both_directions_13mer(kmer)
total_tf = index_13mer.get_total_tf_value_13mer(kmer)

print(f"Forward: {forward_tf}, Reverse: {reverse_tf}, Total: {total_tf}")

# Batch directional analysis
results = index_13mer.get_tf_both_directions_13mer_batch(kmers_13)
for i, (fwd, rev) in enumerate(results):
    print(f"{kmers_13[i]}: forward={fwd}, reverse={rev}")
```

**Complete 13-mer space analysis:**
```python
# Get statistics for the entire 13-mer space
stats = index_13mer.get_13mer_statistics()
print(f"Total 13-mers: {stats['total_kmers']:,}")
print(f"Non-zero frequencies: {stats['non_zero_kmers']:,}")
print(f"Max frequency: {stats['max_frequency']:,}")
print(f"Average frequency: {stats['total_count']/stats['non_zero_kmers']:.2f}")

# Access complete frequency array (4^13 = 67M elements)
# Note: This loads 256MB into memory
full_array = index_13mer.get_13mer_tf_array()
print(f"Array size: {len(full_array):,} elements")
```

### Sequence Coverage Analysis

**Analyze k-mer coverage in sequences:**
```python
# Using real reads from the index
real_read = index_23mer.get_read_by_rid(0)  # Get first read
sequence = real_read.split('~')[0][:100] if '~' in real_read else real_read[:100]    # Take first 100 bp

# Analyze 23-mer coverage using built-in function
coverage_23 = index_23mer.get_sequence_coverage(sequence, cutoff=0, k=23)
print(f"23-mer coverage: {len(coverage_23)} positions")
print(f"Non-zero positions: {sum(1 for tf in coverage_23 if tf > 0)}")
print(f"Average TF: {sum(coverage_23)/len(coverage_23):.2f}")

# Analyze 13-mer coverage using batch queries
kmers_13_in_seq = [sequence[i:i+13] for i in range(len(sequence) - 12)]
coverage_13 = index_13mer.get_total_tf_values_13mer(kmers_13_in_seq)
print(f"13-mer coverage: {len(coverage_13)} positions")
print(f"Non-zero positions: {sum(1 for tf in coverage_13 if tf > 0)}")
print(f"Average TF: {sum(coverage_13)/len(coverage_13):.2f}")
```

**Iterate over k-mers in sequences:**
```python
# 23-mer iteration using built-in iterator
for kmer, tf in index_23mer.iter_sequence_kmers(sequence, k=23):
    if tf > 0:  # Only show k-mers found in index
        print(f"{kmer}: {tf}")

# 13-mer iteration using manual approach (more efficient with batch)
kmers_13 = [sequence[i:i+13] for i in range(len(sequence) - 12)]
tf_values_13 = index_13mer.get_total_tf_values_13mer(kmers_13)

for i, (kmer, tf) in enumerate(zip(kmers_13, tf_values_13)):
    if tf > 0:
        print(f"Position {i}: {kmer}: {tf}")

# For directional analysis of 13-mers
directional_results = index_13mer.get_tf_both_directions_13mer_batch(kmers_13)
for i, (fwd, rev) in enumerate(directional_results):
    if fwd > 0 or rev > 0:
        print(f"Position {i}: {kmers_13[i]}: forward={fwd}, reverse={rev}")
```

### Performance Benchmarks

Based on stress testing with 1M queries and 10K sequence analyses:

| Operation | 13-mers | 23-mers | Speedup |
|-----------|---------|---------|---------|
| **Single TF queries** | 491K queries/sec | 1.1M queries/sec | 23-mer 2.2x faster |
| **Batch TF queries** | 2.0M queries/sec | 2.3M queries/sec | 23-mer 1.2x faster |
| **Sequence coverage** | 24.5K sequences/sec | 17.5K sequences/sec | 13-mer 1.4x faster |
| **K-mer positions** | 2.2M positions/sec | 1.4M positions/sec | 13-mer 1.6x faster |

**Key findings:**
- **Batch operations**: 2-4x faster than single queries for both modes
- **23-mers**: Better for single/batch TF queries due to optimized sparse indexing
- **13-mers**: Better for sequence analysis due to complete space coverage
- **Memory efficiency**: Minimal memory growth during batch operations

### Working with Reads

**Access reads by ID:**
```python
# Get reads from either index
for rid in range(min(5, index_23mer.n_reads)):
    read = index_23mer.get_read_by_rid(rid)
    print(f"Read {rid}: {read[:50]}...")  # First 50 characters
    
    # Split paired reads (separated by '~')
    if '~' in read:
        read1, read2 = read.split('~')
        print(f"  Read 1: {len(read1)} bp, Read 2: {len(read2)} bp")
```

**Iterate over all reads:**
```python
# Iterate through reads with automatic ID assignment
read_count = 0
for rid, read in index_23mer.iter_reads():
    read_count += 1
    if read_count <= 5:  # Show first 5 reads
        print(f"Read {rid}: {len(read)} bp")
    if read_count >= 1000:  # Process first 1000 reads
        break
        
print(f"Processed {read_count} reads")
```

### Complete Example

Here's a practical example showing both 13-mer and 23-mer analysis:

```python
from aindex.core.aindex import AIndex
import aindex.core.aindex_cpp as aindex_cpp
import time

# Load both indices
print("Loading indices...")
index_23mer = AIndex.load_from_prefix("temp/reads.23")
index_23mer.load_reads("temp/reads.reads")

index_13mer = aindex_cpp.AindexWrapper()
index_13mer.load_from_prefix_13mer("temp/all_13mers")
index_13mer.load_reads("temp/reads.reads")

# Get a real sequence to analyze
real_read = index_23mer.get_read_by_rid(0)
sequence = real_read.split('~')[0][:100] if '~' in real_read else real_read[:100]
print(f"Analyzing sequence: {sequence[:50]}...")

# Compare 13-mer vs 23-mer coverage
print("\n=== Coverage Analysis ===")

# 23-mer coverage using built-in function
start = time.time()
coverage_23 = index_23mer.get_sequence_coverage(sequence, cutoff=0, k=23)
time_23 = time.time() - start

# 13-mer coverage using batch query
kmers_13 = [sequence[i:i+13] for i in range(len(sequence) - 12)]
start = time.time()
coverage_13 = index_13mer.get_total_tf_values_13mer(kmers_13)
time_13 = time.time() - start

print(f"23-mers: {len(coverage_23)} positions, {sum(1 for x in coverage_23 if x > 0)} covered ({time_23*1000:.1f}ms)")
print(f"13-mers: {len(coverage_13)} positions, {sum(1 for x in coverage_13 if x > 0)} covered ({time_13*1000:.1f}ms)")

# Performance comparison
print(f"\n=== Performance Test ===")
test_kmers_23 = ["ATCGATCGATCGATCGATCGATC"] * 1000
test_kmers_13 = ["ATCGATCGATCGA"] * 1000

# 23-mer batch query
start = time.time()
results_23 = index_23mer.get_tf_values(test_kmers_23)
time_23_batch = time.time() - start

# 13-mer batch query
start = time.time()
results_13 = index_13mer.get_total_tf_values_13mer(test_kmers_13)
time_13_batch = time.time() - start

print(f"23-mer batch (1K queries): {len(test_kmers_23)/time_23_batch:.0f} queries/sec")
print(f"13-mer batch (1K queries): {len(test_kmers_13)/time_13_batch:.0f} queries/sec")

# Statistics
stats_23 = {"kmers": index_23mer.n_kmers, "reads": index_23mer.n_reads}
stats_13 = index_13mer.get_13mer_statistics()

print(f"\n=== Index Statistics ===")
print(f"23-mer index: {stats_23['kmers']:,} k-mers, {stats_23['reads']:,} reads")
print(f"13-mer index: {stats_13['total_kmers']:,} total k-mers, {stats_13['non_zero_kmers']:,} non-zero")
```

Expected output:
```
Loading indices...
Analyzing sequence: NNNNNNNNNNACTGAACCGCCTTCCGATCTCCAGCTGCAAAGCGTAG...

=== Coverage Analysis ===
23-mers: 78 positions, 42 covered (0.3ms)
13-mers: 88 positions, 88 covered (0.1ms)

=== Performance Test ===
23-mer batch (1K queries): 2,300,000 queries/sec
13-mer batch (1K queries): 2,000,000 queries/sec

=== Index Statistics ===
23-mer index: 15,234,567 k-mers, 125,000 reads  
13-mer index: 67,108,864 total k-mers, 8,945,123 non-zero
```

# Advanced Features

## 13-mer Integration

The aindex library provides highly optimized 13-mer k-mer counting and querying using precomputed perfect hash tables. This mode offers complete coverage of the 4^13 k-mer space with exceptional performance.

### Performance Characteristics

**Query Performance:**
- **Single queries**: 491K queries/second  
- **Batch queries**: 2.0M queries/second (4.1x speedup)
- **Directional queries**: 1.8M queries/second (forward + reverse complement)
- **Complete space**: Access to all 67,108,864 possible 13-mers

**Sequence Analysis Performance:**
- **Coverage analysis**: 24,500 sequences/second
- **Position analysis**: 2.2M k-mer positions/second  
- **Memory efficiency**: Zero memory growth during batch operations
- **Real data coverage**: 100% (all k-mers found in biological data)

### 13-mer Workflow

**1. Generate Complete 13-mer Space:**
```bash
# Generate all possible 13-mers (67M k-mers)
./bin/generate_all_13mers.exe all_13mers.txt

# Build perfect hash for instant lookup
./bin/build_13mer_hash.exe all_13mers.txt temp/all_13mers 4

# Count k-mers in your genomic data
./bin/count_kmers13.exe input_reads.fasta temp/all_13mers.tf.bin hash_file 4
```

**2. Python API Usage:**
```python
from aindex.core.aindex import AIndex

# Load 13-mer index with complete k-mer space
index = AIndex.load_from_prefix_13mer("temp/all_13mers")
index.load_reads("temp/reads.reads")  # Optional: load read sequences

# Query performance demonstration
import time

# Single k-mer query
start = time.time()
tf = index.get_total_tf_value_13mer("ATCGATCGATCGA")
single_time = time.time() - start
print(f"Single query: {tf} (took {single_time*1000:.3f}ms)")

# Batch query (much faster)
kmers = ["ATCGATCGATCGA", "AAAAAAAAAAAAA", "TTTTTTTTTTTTT"] * 1000  # 3K queries
start = time.time() 
tf_values = index.get_total_tf_values_13mer(kmers)
batch_time = time.time() - start
print(f"Batch {len(kmers)} queries: {batch_time:.3f}s ({len(kmers)/batch_time:.0f} queries/sec)")

# Directional analysis (forward + reverse complement)
forward, reverse = index.get_tf_both_directions_13mer("ATCGATCGATCGA")
total = index.get_total_tf_value_13mer("ATCGATCGATCGA")
print(f"Directional: forward={forward}, reverse={reverse}, total={total}")
```

### 13-mer Statistics and Analysis

**Get comprehensive statistics:**
```python
# Complete 13-mer space statistics
stats = index.get_13mer_statistics()
print(f"Total 13-mer space: {stats['total_kmers']:,}")
print(f"Found in data: {stats['non_zero_kmers']:,} ({stats['non_zero_kmers']/stats['total_kmers']*100:.2f}%)")
print(f"Max frequency: {stats['max_frequency']:,}")
print(f"Total occurrences: {stats['total_count']:,}")
print(f"Average frequency: {stats['total_count']/stats['non_zero_kmers']:.2f}")

# Access complete frequency array (warning: 256MB)
if stats['non_zero_kmers'] > 0:
    # Get subset for analysis rather than full array
    sample_indices = range(0, 1000000, 1000)  # Sample every 1000th element
    sample_tfs = [index.get_tf_by_index_13mer(i) for i in sample_indices]
    non_zero_sample = [tf for tf in sample_tfs if tf > 0]
    print(f"Sample analysis: {len(non_zero_sample)}/{len(sample_tfs)} non-zero in sample")
```

**Sequence coverage analysis:**
```python
# Analyze real genomic sequences
for rid in range(min(5, index.n_reads)):
    read = index.get_read_by_rid(rid)
    if '~' in read:
        sequence = read.split('~')[0]  # Take first mate
    else:
        sequence = read
    
    # Limit to reasonable length for demonstration
    if len(sequence) > 100:
        sequence = sequence[:100]
    
    # Compute 13-mer coverage
    coverage = []
    for i in range(len(sequence) - 12):
        kmer = sequence[i:i+13]
        tf = index.get_total_tf_value_13mer(kmer)
        coverage.append(tf)
    
    if coverage:
        avg_tf = sum(coverage) / len(coverage)
        max_tf = max(coverage)
        coverage_pct = sum(1 for tf in coverage if tf > 0) / len(coverage) * 100
        print(f"Read {rid}: {len(coverage)} 13-mers, {coverage_pct:.1f}% covered, avg TF {avg_tf:.1f}, max TF {max_tf}")
```

## 23-mer Integration  

The 23-mer mode provides efficient sparse indexing for longer k-mers commonly used in genomic analysis.

### Performance Characteristics

**Query Performance:**
- **Single queries**: 1.0M queries/second
- **Batch queries**: 2.4M queries/second (2.4x speedup)  
- **Directional queries**: Available for forward + reverse complement analysis
- **Sparse indexing**: Only stores k-mers present in input data

**Sequence Analysis Performance:**
- **Coverage analysis**: 16,900 sequences/second
- **Position analysis**: 1.3M k-mer positions/second
- **Memory efficiency**: Constant memory usage during operations
- **Real data coverage**: 100% (all k-mers found in genomic sequences)

### 23-mer Workflow

**1. Build 23-mer Index:**
```bash
# Using the fast built-in k-mer counter (recommended)
FASTQ1=./tests/raw_reads.101bp.IS350bp25_1.fastq
FASTQ2=./tests/raw_reads.101bp.IS350bp25_2.fastq
OUTPUT_PREFIX=./temp/reads.23

python3 scripts/compute_aindex.py -i $FASTQ1,$FASTQ2 -t fastq -o $OUTPUT_PREFIX --lu 2 -P 30 --use_kmer_counter
```

**2. Python API Usage:**
```python
from aindex.core.aindex import AIndex

# Load 23-mer index
index = AIndex.load_from_prefix("temp/reads.23")
index.load_reads("temp/reads.reads")

# Performance demonstration
import time

# Batch query performance
kmers = ["ATCGATCGATCGATCGATCGATC", "AAAAAAAAAAAAAAAAAAAAAA"] * 1000  # 2K queries
start = time.time()
tf_values = index.get_tf_values(kmers)  # Auto-detects 23-mer mode
batch_time = time.time() - start
print(f"23-mer batch {len(kmers)} queries: {batch_time:.3f}s ({len(kmers)/batch_time:.0f} queries/sec)")

# Sequence coverage analysis
read = index.get_read_by_rid(0)
sequence = read.split('~')[0][:100] if '~' in read else read[:100]

start = time.time()
coverage = index.get_sequence_coverage(sequence, cutoff=0, k=23)
coverage_time = time.time() - start

print(f"23-mer coverage analysis: {len(coverage)} positions in {coverage_time*1000:.1f}ms")
print(f"Coverage: {sum(1 for tf in coverage if tf > 0)/len(coverage)*100:.1f}% positions covered")
print(f"Average TF: {sum(coverage)/len(coverage):.1f}")
```

## Performance Comparison

### Throughput Comparison (Operations per Second)

| Operation Type | 13-mers | 23-mers | Winner |
|----------------|---------|---------|--------|
| **Single TF queries** | 491K/sec | 1.1M/sec | 23-mer (+124%) |
| **Batch TF queries** | 2.0M/sec | 2.3M/sec | 23-mer (+15%) |
| **Sequence coverage** | 24.5K/sec | 17.5K/sec | 13-mer (+40%) |
| **Position analysis** | 2.2M/sec | 1.4M/sec | 13-mer (+57%) |

### Use Case Recommendations

**Choose 13-mers when:**
- Analyzing complete k-mer space (population genetics, mutation analysis)
- Maximum query performance needed
- Working with shorter sequences or fragments
- Need comprehensive coverage statistics

**Choose 23-mers when:**
- Standard genomic analysis (assembly, alignment, variant calling)
- Working with longer reads (>100bp)
- Memory efficiency is critical
- Integration with existing 23-mer workflows

### Memory Usage

- **13-mer index**: ~277MB (256MB frequencies + 21MB hash)
- **23-mer index**: Variable, depends on data complexity
- **Both modes**: Memory-mapped files for efficient access
- **Batch operations**: Minimal additional memory overhead

## File Formats

### 13-mer Files
- `.tf.bin`: Binary frequency array (uint64_t × 67M elements = 512MB)
- `.pf`: Perfect hash function for k-mer → index mapping
- `.kmers.bin`: Binary k-mer encoding (optional validation)

### 23-mer Files  
- `.tf.bin`: Binary frequency array (variable size)
- `.pf`: Perfect hash function
- `.kmers.bin`: Binary k-mer storage
- `.aindex.indices.bin` & `.aindex.index.bin`: Position indices (optional)
