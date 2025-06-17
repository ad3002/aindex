import sys
import os
import logging
import numpy as np
import time
import subprocess

# Path setup - add project root directory to Python path
PATH_TO_AINDEX_FOLDER = '.'
sys.path.insert(0, PATH_TO_AINDEX_FOLDER)

# Logging setup
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Import aindex module - direct import of C++ module for 13-mers
try:
    import aindex.core.aindex_cpp as aindex_cpp
    print("✓ aindex_cpp module imported successfully")
    
    # Try to import Python wrapper as well
    try:
        from aindex.core.aindex import AIndex
        print("✓ AIndex module also imported successfully")
        HAS_PYTHON_WRAPPER = True
    except ImportError as e:
        print(f"⚠️ Python wrapper unavailable: {e}")
        HAS_PYTHON_WRAPPER = False
        AIndex = None
        
except ImportError as e:
    print(f"✗ Import error for aindex_cpp: {e}")
    print("Check that the aindex package is installed and built correctly")
    print("Run in terminal:")
    print(f"  cd {PATH_TO_AINDEX_FOLDER}")
    print("  make pybind11")
    aindex_cpp = None
    AIndex = None
    HAS_PYTHON_WRAPPER = False

print("Working directory:", os.getcwd())
print("Python version:", sys.version)
print("Path to aindex:", PATH_TO_AINDEX_FOLDER)

assert aindex_cpp is not None, "aindex_cpp module not imported, tests cannot be run"

# Check for required files before loading
prefix = os.path.join(PATH_TO_AINDEX_FOLDER, "temp/all_13mers")
required_files = [
    f"{prefix}.pf",
    f"{prefix}.index.bin", 
    f"{prefix}.kmers.bin"
]

# Check for reads files (path without k-mer size)
base_prefix = os.path.join(PATH_TO_AINDEX_FOLDER, "temp/reads")
reads_file = f"{base_prefix}.reads"
reads_index_file = f"{base_prefix}.ridx"  # Correct path to reads index
kmers_file = f"{prefix}.kmers"
has_reads = os.path.exists(reads_file) and os.path.exists(reads_index_file)

# Function to count lines in a file
def count_lines(file_path):
    try:
        result = subprocess.run(["wc", "-l", file_path], capture_output=True, text=True)
        if result.returncode == 0:
            return int(result.stdout.split()[0])
        return None
    except Exception as e:
        print(f"Error counting lines in {file_path}: {e}")
        return None

def parse_kmers_analysis_file(file_path):
    """
    Parses kmers_analysis.trues file and returns dict: kmer -> frequency
    
    File format:
    kmer    frequency   [read_id,position,direction] [read_id,position,direction] ...
    """
    kmer_tf = {}
    
    if not os.path.exists(file_path):
        print(f"File {file_path} not found, skipping k-mer check")
        return kmer_tf
        
    try:
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    kmer = parts[0]
                    frequency = int(parts[1])
                    kmer_tf[kmer] = frequency
        print(f"Loaded {len(kmer_tf)} kmers from file {file_path}")
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        
    return kmer_tf

# Get expected number of reads and kmers
expected_reads_count = count_lines(reads_file) if os.path.exists(reads_file) else None
expected_kmers_count = count_lines(kmers_file) if os.path.exists(kmers_file) else None

print("Checking for required files:")
all_files_exist = True
for file_path in required_files:
    exists = os.path.exists(file_path)
    print(f"  {file_path}: {'✓' if exists else '✗'}")
    if not exists:
        all_files_exist = False

print(f"  {reads_file}: {'✓' if os.path.exists(reads_file) else '✗'}")
print(f"  {reads_index_file}: {'✓' if os.path.exists(reads_index_file) else '✗'}")
print(f"  {kmers_file}: {'✓' if os.path.exists(kmers_file) else '✗'}")

print("\nExpected counts:")
if expected_reads_count:
    print(f"  Number of reads (wc -l {reads_file}): {expected_reads_count:,}")
if expected_kmers_count:
    print(f"  Number of k-mers (wc -l {kmers_file}): {expected_kmers_count:,}")

if not all_files_exist:
    print("\n⚠️ Not all required files found. Check that the index was created correctly.")
    print(f"Expected files should be in directory: {os.path.dirname(prefix)}")
else:
    print("\n✓ All required files found, starting loading...")

# Load index only if aindex_cpp is available
if aindex_cpp is not None:
    try:
        start_time = time.time()
        print(f"Loading 13-mer index from prefix: {prefix}")
        
        # Create AindexWrapper instance and load data
        index = aindex_cpp.AindexWrapper()
        
        # Load 13-mer index
        index.load_from_prefix_13mer(prefix)
        
        # Load reads if available
        if has_reads:
            print("Loading reads into 13-mer index...")
            index.load_reads(reads_file)
            print(f"✓ Reads loaded: {index.n_reads:,} reads")
        else:
            print("⚠️ Reads files not found, loading only TF data")
        
        load_time = time.time() - start_time
        print(f"✓ 13-mer index loaded successfully in {load_time:.2f} seconds")
        
        # Basic info about 13-mer index
        stats = index.get_13mer_statistics()
        print(f"Total 13-mers: {stats['total_kmers']:,}")
        print(f"13-mers with non-zero TF: {stats['non_zero_kmers']:,}")
        print(f"Max frequency: {stats['max_frequency']:,}")
        print(f"Total TF count: {stats['total_count']:,}")
        print(f"Number of reads: {index.n_reads:,}")
        print(f"Reads file size: {index.reads_size:,} bytes")
        
        print(f"Index loaded for 13-mers: True")
        
        # Check that reads were loaded
        if index.n_reads > 0:
            print(f"✓ Reads loaded successfully ({index.n_reads:,} reads)")
            
            # Demonstrate access to reads
            try:
                # Get first read for demonstration
                first_read = index.get_read_by_rid(0)
                if first_read:
                    print(f"First read (first 50 chars): {first_read[:50]}...")
            except Exception as e:
                print(f"⚠️ Error accessing reads: {e}")
        else:
            print("⚠️ Reads not loaded or count is 0")
        
        # Remove expected vs actual count check (using C++ API only)
            
        # TF file statistics for 13-mers
        print("\n=== TF file statistics for 13-mers ===")
        try:
            # Get statistics via new API
            stats = index.get_13mer_statistics()
            total_kmers = stats['total_kmers']
            non_zero_count = stats['non_zero_kmers']
            max_tf = stats['max_frequency']
            total_count = stats['total_count']
            avg_tf = total_count / non_zero_count if non_zero_count > 0 else 0
            
            print(f"Total 13-mers: {total_kmers:,}")
            print(f"13-mers with non-zero TF: {non_zero_count:,} ({non_zero_count*100/total_kmers:.2f}%)")
            print(f"Max TF: {max_tf:,}")
            print(f"Total TF count: {total_count:,}")
            print(f"Average TF (non-zero): {avg_tf:.2f}")
        except Exception as e:
            print(f"Error getting TF statistics: {e}")
            
        # Remove reads check, using only TF API for 13-mers
        
        # Check k-mers and their frequencies against expected values
        kmers_analysis_file = os.path.join(PATH_TO_AINDEX_FOLDER, "./temp/all_13mers.true.kmers")
        if os.path.exists(kmers_analysis_file):
            print("\n=== Checking 13-mers against expected values ===")
            expected_kmers = parse_kmers_analysis_file(kmers_analysis_file)
            
            if expected_kmers:
                print(f"Loaded {len(expected_kmers)} expected 13-mers from {kmers_analysis_file}")
                
                # Test first few k-mers for demonstration
                sample_size = min(10, len(expected_kmers))
                sample_kmers = list(expected_kmers.keys())[:sample_size]
                
                print(f"\nChecking {sample_size} random 13-mers:")
                print("Kmer\t\t\tExpected\tForward\tReverse\tTotal\tStatus")
                print("-" * 80)
                
                all_matched = True
                mismatches = 0
                
                for kmer in sample_kmers:
                    expected_tf = expected_kmers[kmer]
                    
                    # Use new functions to get TF in both directions
                    forward_tf, reverse_tf = index.get_tf_both_directions_13mer(kmer)
                    total_tf = index.get_total_tf_value_13mer(kmer)
                    
                    matches = expected_tf == total_tf
                    all_matched = all_matched and matches
                    if not matches:
                        mismatches += 1
                        
                    status = "✓" if matches else "✗"
                    print(f"{status} {kmer}\t{expected_tf}\t\t{forward_tf}\t{reverse_tf}\t{total_tf}\t{status}")
                
                # Check all k-mers
                print("\nChecking all 13-mers:")
                total_mismatches = 0
                zero_tf_count = 0
                
                for kmer, expected_tf in expected_kmers.items():
                    total_tf = index.get_total_tf_value_13mer(kmer)
                    
                    if total_tf == 0:
                        zero_tf_count += 1
                    
                    if expected_tf != total_tf:
                        total_mismatches += 1
                        
                # Output results
                match_percentage = 100 * (len(expected_kmers) - total_mismatches) / len(expected_kmers)
                print(f"Matched: {match_percentage:.2f}% 13-mers")
                print(f"Mismatches: {total_mismatches:,} out of {len(expected_kmers):,} 13-mers")
                print(f"13-mers with zero frequency: {zero_tf_count:,}")
                
                if total_mismatches > 0:
                    print("\n⚠️ Mismatches found between expected and actual 13-mer frequencies")
                else:
                    print("\n✓ All 13-mer frequencies match expected values")
            else:
                print("Failed to load expected 13-mer values for checking")
        else:
            print(f"\n⚠️ File {kmers_analysis_file} not found, skipping 13-mer check")
            print("  To create the file run:\n  python tests/analyze_kmers.py --input-file temp/reads.reads -o kmers_analysis.trues")
    except Exception as e:
        print(f"✗ Error loading index: {e}")
        import traceback
        traceback.print_exc()
        index = None
else:
    print("⚠️ aindex_cpp not imported, skipping index loading")
    index = None

# Get detailed info about the index
if index is None:
    print("⚠️ Index not loaded, skipping testing")
else:
    print("=== 13-mer index information ===")
    try:
        stats = index.get_13mer_statistics()
        print(f"Total 13-mers: {stats['total_kmers']:,}")
        print(f"13-mers with non-zero TF: {stats['non_zero_kmers']:,}")
        print(f"Max frequency: {stats['max_frequency']:,}")
        print(f"Total TF count: {stats['total_count']:,}")
    except Exception as e:
        print(f"Error getting statistics: {e}")

    print("\n=== Additional statistics ===")
    print(f"13-mer mode: active")

    # Test access to TF for one 13-mer
    test_kmer = "AAAAACCCCCGGG"  # 13-mer for testing
    print(f"\nTF access test:")
    print(f"  13-mer: '{test_kmer}'")
    try:
        forward_tf, reverse_tf = index.get_tf_both_directions_13mer(test_kmer)
        total_tf = index.get_total_tf_value_13mer(test_kmer)
        rc_kmer = index.get_reverse_complement_13mer(test_kmer)
        
        print(f"  Forward TF: {forward_tf}")
        print(f"  Reverse TF: {reverse_tf} (RC: {rc_kmer})")
        print(f"  Total TF: {total_tf}")
    except Exception as e:
        print(f"  Error: {e}")

# Test getting term frequency for 13-mers
if index is None:
    print("⚠️ Index not loaded, skipping TF testing")
else:
    print("=== Test getting TF for 13-mers ===")

    # Create several test 13-mers
    test_kmers = [
        "AAAAAAAAAAAAA",   # only A
        "TTTTTTTTTTTTT",   # only T
        "GGGGGGGGGGGGG",   # only G
        "CCCCCCCCCCCCC",   # only C
        "ATCGATCGATCGA",   # random sequence
        "AAAAAGAGTTAAT",   # known k-mer from our tests
        "AGTAGTAGTAGTA"    # another known k-mer
    ]

    print("Testing individual 13-mers:")
    print("Kmer\t\t\tForward\tReverse\tTotal\tRC")
    print("-" * 80)
    
    for kmer in test_kmers:
        try:
            forward_tf, reverse_tf = index.get_tf_both_directions_13mer(kmer)
            total_tf = index.get_total_tf_value_13mer(kmer)
            rc_kmer = index.get_reverse_complement_13mer(kmer)
            print(f"{kmer}\t{forward_tf}\t{reverse_tf}\t{total_tf}\t{rc_kmer}")
        except Exception as e:
            print(f"{kmer}\tError: {e}")

    print("\n=== Test batch getting TF for 13-mers ===")

    # Test batch processing
    batch_kmers = test_kmers[:5]  # take first 5
    try:
        start_time = time.time()
        
        # Batch get TF in both directions
        batch_results = index.get_tf_both_directions_13mer_batch(batch_kmers)
        
        # Batch get total TF
        total_tfs = index.get_total_tf_values_13mer(batch_kmers)
        
        batch_time = time.time() - start_time
        
        print(f"Batch processing {len(batch_kmers)} 13-mers took {batch_time:.4f} seconds")
        print("Results:")
        print("Kmer\t\t\tForward\tReverse\tTotal")
        print("-" * 60)
        
        for i, kmer in enumerate(batch_kmers):
            forward_tf, reverse_tf = batch_results[i]
            total_tf = total_tfs[i]
            print(f"{kmer}\t{forward_tf}\t{reverse_tf}\t{total_tf}")
            
    except Exception as e:
        print(f"Batch processing error: {e}")

    print("\n=== Performance comparison for 13-mers ===")

    # Compare single calls vs batch
    n_test = 100
    random_kmers = []
    for i in range(n_test):
        kmer = ''.join(np.random.choice(['A', 'T', 'G', 'C'], 13))
        random_kmers.append(kmer)

    # Single calls
    start_time = time.time()
    single_results = []
    for kmer in random_kmers:
        try:
            total_tf = index.get_total_tf_value_13mer(kmer)
            single_results.append(total_tf)
        except:
            single_results.append(0)
    single_time = time.time() - start_time

    # Batch call
    start_time = time.time()
    try:
        batch_results = index.get_total_tf_values_13mer(random_kmers)
    except:
        batch_results = [0] * len(random_kmers)
    batch_time = time.time() - start_time

    print(f"Single calls ({n_test} 13-mers): {single_time:.4f} sec")
    print(f"Batch call ({n_test} 13-mers): {batch_time:.4f} sec")
    print(f"Speedup: {single_time/batch_time:.1f}x")

# For 13-mer API, positions and reads are not available, main testing ends here
print("\n=== Test analysis of random 13-mers ===")

# ... (rest of the code continues in English, following the same translation pattern)
