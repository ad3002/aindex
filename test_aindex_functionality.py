import sys
import os
import logging
import numpy as np
import time
import subprocess

# Path setup - add the project root directory to Python path
PATH_TO_AINDEX_FOLDER = '..'
sys.path.insert(0, PATH_TO_AINDEX_FOLDER)

# Logging setup
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Import aindex module
try:
    from aindex.core.aindex import AIndex
    print("✓ AIndex module imported successfully")
except ImportError as e:
    print(f"✗ Import error for AIndex: {e}")
    print("Check that the aindex package is installed and built correctly")
    print("Run in terminal:")
    print(f"  cd {PATH_TO_AINDEX_FOLDER}")
    print("  python setup.py build_ext --inplace")
    AIndex = None

print("Working directory:", os.getcwd())
print("Python version:", sys.version)
print("Path to aindex:", PATH_TO_AINDEX_FOLDER)

assert AIndex is not None, "AIndex module not imported, tests cannot be run"

# Check for required files before loading
prefix = os.path.join(PATH_TO_AINDEX_FOLDER, "temp/reads.23")
required_files = [
    f"{prefix}.pf",
    f"{prefix}.tf.bin", 
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
    Parses kmers_analysis.trues file and returns a dict: kmer -> frequency
    
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
        print(f"Loaded {len(kmer_tf)} kmers from {file_path}")
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        
    return kmer_tf

# Get expected number of reads and k-mers
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

# Load index only if AIndex is available
if AIndex is not None:
    try:
        start_time = time.time()
        print(f"Loading index from prefix: {prefix}")
        
        # Load index with auto-detection mode and reads (if available)
        if has_reads:
            print("Loading index with reads...")
            # Pass reads_file directly, since load_reads=True looks for reads_file with same prefix
            index = AIndex.load_from_prefix(prefix)
            index.load_reads(reads_file)
        else:
            print("Loading index without reads (reads files not found)...")
            index = AIndex.load_from_prefix(prefix)
        
        load_time = time.time() - start_time
        print(f"✓ Index loaded successfully in {load_time:.2f} seconds")
        
        # Basic info
        print(f"Number of k-mers: {index.n_kmers:,}")
        print(f"Number of reads: {index.n_reads:,}")
        print(f"AIndex loaded: {index.aindex_loaded}")
        
        # Compare expected and actual counts
        if expected_kmers_count:
            print(f"K-mer check: {'✓' if index.n_kmers == expected_kmers_count else '✗'} " +
                 f"(loaded {index.n_kmers:,}, expected {expected_kmers_count:,})")
            
        if has_reads and expected_reads_count:
            print(f"Reads check: {'✓' if index.n_reads == expected_reads_count else '✗'} " +
                 f"(loaded {index.n_reads:,}, expected {expected_reads_count:,})")
        
        # Check that reads were loaded
        if has_reads and index.n_reads > 0:
            print(f"✓ Reads loaded successfully ({index.n_reads:,} reads)")
            
            # Demonstrate access to reads if loaded
            try:
                # Get first read for demonstration
                first_read = index.get_read_by_rid(0)
                if first_read:
                    print(f"First read (first 50 bases): {first_read[:50]}...")
            except Exception as e:
                print(f"⚠️ Error accessing reads: {e}")
        else:
            print("⚠️ Reads not loaded or count is 0")
        
        # Check k-mers and their frequencies against expected values
        kmers_analysis_file = os.path.join(PATH_TO_AINDEX_FOLDER, "kmers_analysis.trues")
        if os.path.exists(kmers_analysis_file):
            print("\n=== Checking k-mers against expected values ===")
            expected_kmers = parse_kmers_analysis_file(kmers_analysis_file)
            
            if expected_kmers:
                print(f"Loaded {len(expected_kmers)} expected k-mers from {kmers_analysis_file}")
                
                # Test first few k-mers for demonstration
                sample_size = min(10, len(expected_kmers))
                sample_kmers = list(expected_kmers.keys())[:sample_size]
                
                print(f"\nChecking {sample_size} random k-mers:")
                all_matched = True
                mismatches = 0
                
                for kmer in sample_kmers:
                    expected_tf = expected_kmers[kmer]
                    actual_tf = index[kmer]
                    
                    matches = expected_tf == actual_tf
                    all_matched = all_matched and matches
                    if not matches:
                        mismatches += 1
                        
                    status = "✓" if matches else "✗"
                    print(f"{status} {kmer}: expected {expected_tf}, got {actual_tf}")
                
                # Check all k-mers
                print("\nChecking all k-mers:")
                total_mismatches = 0
                zero_tf_count = 0
                
                for kmer, expected_tf in expected_kmers.items():
                    actual_tf = index[kmer]
                    
                    if actual_tf == 0:
                        zero_tf_count += 1
                    
                    if expected_tf != actual_tf:
                        total_mismatches += 1
                        
                # Output results
                match_percentage = 100 * (len(expected_kmers) - total_mismatches) / len(expected_kmers)
                print(f"Matched: {match_percentage:.2f}% k-mers")
                print(f"Mismatches: {total_mismatches:,} out of {len(expected_kmers):,} k-mers")
                print(f"K-mers with zero frequency: {zero_tf_count:,}")
                
                if total_mismatches > 0:
                    print("\n⚠️ Mismatches found between expected and actual k-mer frequencies")
                else:
                    print("\n✓ All k-mer frequencies match expected values")
            else:
                print("Could not load expected k-mer values for checking")
        else:
            print(f"\n⚠️ File {kmers_analysis_file} not found, skipping k-mer check")
            print("  To create the file, run:\n  python tests/analyze_kmers.py --input-file temp/reads.reads -o kmers_analysis.trues")
    except Exception as e:
        print(f"✗ Error loading index: {e}")
        import traceback
        traceback.print_exc()
        index = None
else:
    print("⚠️ AIndex not imported, skipping index loading")
    index = None

# Get detailed index info
if index is None:
    print("⚠️ Index not loaded, skipping testing")
else:
    print("=== Index information ===")
    try:
        info = index.get_index_info()
        print(info)
    except Exception as e:
        print(f"Error getting info: {e}")

    print("\n=== Additional statistics ===")
    print(f"Hash size: {index.get_hash_size():,}")
    print(f"Reads data size: {index.get_reads_size():,}")

    # Check main properties
    print(f"\nIndex properties:")
    print(f"  n_kmers: {index.n_kmers:,}")
    print(f"  n_reads: {index.n_reads:,}")
    print(f"  aindex_loaded: {index.aindex_loaded}")

    # Check magic methods
    print(f"\nMagic methods:")
    print(f"  len(index): {len(index):,}")

    # Test index as a dictionary
    test_kmer = "AAAAAATGAATAACGTTATTATT"
    print(f"\nTest index access:")
    print(f"  index['{test_kmer}']: ", end="")
    try:
        tf_value = index[test_kmer]
        print(f"{tf_value}")
    except Exception as e:
        print(f"Error: {e}")

# Test getting term frequency for a single k-mer
if index is None:
    print("⚠️ Index not loaded, skipping TF testing")
else:
    print("=== Test getting TF for a single k-mer ===")

    # Create several test 23-mers
    test_kmers = [
        "ATCGATCGATCGATCGATCGATC",  # all same nucleotides
        "AAAAAAAAAAAAAAAAAAAAAA",   # only A
        "TTTTTTTTTTTTTTTTTTTTTT",   # only T
        "GGGGGGGGGGGGGGGGGGGGGGG",  # only G
        "CCCCCCCCCCCCCCCCCCCCCCC",  # only C
        "ATCGATCGATCGATCGATCGATG",  # random sequence
        "NNNNNNNNNNNNNNNNNNNNNNN"   # undefined nucleotides
    ]

    for kmer in test_kmers:
        try:
            tf = index[kmer]
            print(f"  {kmer}: tf = {tf}")
        except Exception as e:
            print(f"  {kmer}: Error = {e}")

    print("\n=== Test getting TF for multiple k-mers at once ===")

    # Test batch processing
    batch_kmers = test_kmers[:5]  # take first 5
    try:
        start_time = time.time()
        tf_values = index.get_tf_values(batch_kmers)
        batch_time = time.time() - start_time
        
        print(f"Batch processing {len(batch_kmers)} k-mers took {batch_time:.4f} seconds")
        print("Results:")
        for kmer, tf in zip(batch_kmers, tf_values):
            print(f"  {kmer}: tf = {tf}")
            
    except Exception as e:
        print(f"Batch processing error: {e}")

    print("\n=== Performance comparison ===")

    # Compare single calls vs batch
    n_test = 100
    random_kmers = []
    for i in range(n_test):
        kmer = ''.join(np.random.choice(['A', 'T', 'G', 'C'], 23))
        random_kmers.append(kmer)

    # Single calls
    start_time = time.time()
    single_results = []
    for kmer in random_kmers:
        try:
            tf = index.get_tf_value(kmer)
            single_results.append(tf)
        except:
            single_results.append(0)
    single_time = time.time() - start_time

    # Batch call
    start_time = time.time()
    try:
        batch_results = index.get_tf_values(random_kmers)
    except:
        batch_results = [0] * len(random_kmers)
    batch_time = time.time() - start_time

    print(f"Single calls ({n_test} k-mers): {single_time:.4f} sec")
    print(f"Batch call ({n_test} k-mers): {batch_time:.4f} sec")
    print(f"Speedup: {single_time/batch_time:.1f}x")

# Test getting positions of k-mers in reads
if index is None:
    print("⚠️ Index not loaded, skipping position testing")
    test_kmers_for_positions = []
else:
    print("=== Test getting positions of k-mers ===")

    if not index.aindex_loaded:
        print("⚠️ AIndex not loaded, positions unavailable")
        test_kmers_for_positions = []
    else:
        # Find k-mers with nonzero TF for testing
        test_kmers_for_positions = []
        
        print("Searching for k-mers with nonzero TF...")
        for i in range(1000):  # try first 1000 random k-mers
            kmer = ''.join(np.random.choice(['A', 'T', 'G', 'C'], 23))
            tf = index.get_tf_value(kmer)
            if tf > 0:
                test_kmers_for_positions.append((kmer, tf))
                if len(test_kmers_for_positions) >= 5:  # find 5 k-mers
                    break
        
        if not test_kmers_for_positions:
            print("No k-mers with nonzero TF found")
        else:
            print(f"Found {len(test_kmers_for_positions)} k-mers with nonzero TF")
            
            for kmer, tf in test_kmers_for_positions:
                print(f"\nTesting k-mer: {kmer} (TF = {tf})")
                
                try:
                    positions = index.get_positions(kmer)
                    print(f"  Positions found: {len(positions)}")
                    
                    if positions:
                        print(f"  First 10 positions: {positions[:10]}")
                        
                        # Check that number of positions matches TF
                        if len(positions) != tf:
                            print(f"  ⚠️ Mismatch: TF={tf}, but positions={len(positions)}")
                        else:
                            print(f"  ✓ TF matches number of positions")
                            
                except Exception as e:
                    print(f"  Error getting positions: {e}")
                    
        # Also test pos() method (alias for get_positions)
        if test_kmers_for_positions:
            test_kmer, test_tf = test_kmers_for_positions[0]
            print(f"\n=== Test pos() method ===")
            print(f"Testing k-mer: {test_kmer}")
            
            try:
                positions1 = index.get_positions(test_kmer)
                positions2 = index.pos(test_kmer)
                
                print(f"get_positions(): {len(positions1)} positions")
                print(f"pos(): {len(positions2)} positions")
                print(f"Results identical: {positions1 == positions2}")
                
            except Exception as e:
                print(f"Error testing pos(): {e}")

# Test getting reads by read id
if index is None:
    print("⚠️ Index not loaded, skipping reads testing")
else:
    print("=== Test getting reads by read id ===")

    if index.n_reads == 0:
        print("⚠️ Reads not loaded")
    else:
        print(f"Total reads: {index.n_reads:,}")
        
        # Test first few reads
        test_rids = [0, 1, 2, min(10, index.n_reads-1), min(100, index.n_reads-1)]
        
        for rid in test_rids:
            if rid >= index.n_reads:
                continue
                
            try:
                read = index.get_read_by_rid(rid)
                if read:
                    print(f"Read {rid}: length = {len(read)}, first 50 nt: {read[:50]}...")
                    
                    # Check for subread delimiters
                    if '~' in read:
                        subreads = read.split('~')
                        print(f"  Contains {len(subreads)} subreads")
                    else:
                        print(f"  Single read without delimiters")
                        
                else:
                    print(f"Read {rid}: empty or not found")
                    
            except Exception as e:
                print(f"Read {rid}: Error = {e}")

        # Test iterating over reads
        print(f"\n=== Test iterating over reads ===")
        try:
            read_count = 0
            total_length = 0
            subread_count = 0
            
            # Limit number of reads for testing
            max_test_reads = min(1000, index.n_reads)
            
            for rid, read in index.iter_reads():
                if rid >= max_test_reads:
                    break
                    
                read_count += 1
                total_length += len(read)
                subread_count += len(read.split('~'))
                
            print(f"Tested reads: {read_count}")
            print(f"Average read length: {total_length/read_count:.1f}")
            print(f"Average subreads per read: {subread_count/read_count:.1f}")
            
        except Exception as e:
            print(f"Error iterating over reads: {e}")

        # Test iterating over subreads
        print(f"\n=== Test iterating over subreads ===")
        try:
            subread_test_count = 0
            total_subread_length = 0
            max_test_subreads = 100
            
            for rid, subread_idx, subread in index.iter_reads_se():
                if subread_test_count >= max_test_subreads:
                    break
                    
                subread_test_count += 1
                total_subread_length += len(subread)
                
                if subread_test_count <= 5:  # show first 5
                    print(f"  Read {rid}, subread {subread_idx}: length = {len(subread)}, seq = {subread[:30]}...")
            
            if subread_test_count > 0:
                print(f"Tested subreads: {subread_test_count}")
                print(f"Average subread length: {total_subread_length/subread_test_count:.1f}")
                
        except Exception as e:
            print(f"Error iterating over subreads: {e}")

# Test getting sequence by coordinates
if index is None:
    print("⚠️ Index not loaded, skipping sequence retrieval testing")
else:
    print("=== Test getting sequence by position ===")

    if index.reads_size == 0:
        print("⚠️ Reads file not loaded")
    else:
        print(f"Reads file size: {index.reads_size:,} bytes")
        
        # Test various coordinates
        test_positions = [
            (0, 50),          # start of file
            (100, 150),       # middle
            (1000, 1100),     # another position
            (index.reads_size-100, index.reads_size-50) if index.reads_size > 100 else (0, 50)  # end of file
        ]
        
        for start, end in test_positions:
            if start >= index.reads_size or end >= index.reads_size:
                continue
                
            print(f"\nPositions {start}-{end}:")
            
            try:
                # Get sequence without reverse complement
                seq = index.get_read(start, end, revcomp=False)
                print(f"  Forward: length = {len(seq)}, seq = {seq[:50]}...")
                
                # Get sequence with reverse complement
                seq_rc = index.get_read(start, end, revcomp=True)
                print(f"  Reverse complement: length = {len(seq_rc)}, seq = {seq_rc[:50]}...")
                
                # Check that lengths are equal
                if len(seq) == len(seq_rc):
                    print(f"  ✓ Lengths match")
                else:
                    print(f"  ✗ Lengths do not match: {len(seq)} vs {len(seq_rc)}")
                    
            except Exception as e:
                print(f"  Error: {e}")

        # Test correspondence between k-mer positions and sequences
        if index.aindex_loaded and 'test_kmers_for_positions' in globals() and test_kmers_for_positions:
            print(f"\n=== Checking correspondence of positions and sequences ===")
            
            test_kmer, test_tf = test_kmers_for_positions[0]
            positions = index.get_positions(test_kmer)
            
            if positions:
                print(f"Testing k-mer: {test_kmer}")
                print(f"Positions found: {len(positions)}")
                
                # Check first few positions
                for i, pos in enumerate(positions[:5]):
                    try:
                        # Get sequence of k-mer length
                        seq = index.get_read(pos, pos + 23, revcomp=False)
                        
                        print(f"  Position {pos}: {seq}")
                        
                        # Check match
                        if seq == test_kmer:
                            print(f"    ✓ Exact match")
                        else:
                            # Check reverse complement
                            seq_rc = index.get_read(pos, pos + 23, revcomp=True)
                            if seq_rc == test_kmer:
                                print(f"    ✓ Reverse complement match")
                            else:
                                print(f"    ✗ Mismatch (expected: {test_kmer})")
                                
                    except Exception as e:
                        print(f"  Position {pos}: Error = {e}")

        # Test getting rid and start by position
        if index.aindex_loaded:
            print(f"\n=== Test getting RID and start by position ===")
            
            test_positions_rid = [100, 1000, 5000, 10000]
            
            for pos in test_positions_rid:
                if pos >= index.reads_size:
                    continue
                    
                try:
                    rid = index.get_rid(pos)
                    start = index.get_start(pos)
                    
                    print(f"Position {pos}: RID = {rid}, start = {start}")
                    
                    if rid > 0:
                        # Check that we can get this read
                        read = index.get_read_by_rid(rid)
                        if read:
                            print(f"  Read {rid}: length = {len(read)}")
                        else:
                            print(f"  Read {rid}: not found")
                            
                except Exception as e:
                    print(f"Position {pos}: Error = {e}")

# Test sequence coverage analysis by k-mers
if index is None:
    print("⚠️ Index not loaded, skipping sequence coverage testing")
else:
    print("=== Test sequence coverage analysis ===")

    # Create test sequences
    test_sequences = [
        
    ]

    # If we have real reads, take one of them
    if index.n_reads > 0:
        try:
            real_read = index.get_read_by_rid(0)
            if real_read and len(real_read) >= 50:
                # Take first subread if present
                if '~' in real_read:
                    real_sequence = real_read.split('~')[0]
                else:
                    real_sequence = real_read
                
                # Limit length for demonstration
                if len(real_sequence) > 100:
                    real_sequence = real_sequence[:100]
                    
                test_sequences.append(real_sequence)
                print(f"Added real sequence of length {len(real_sequence)}")
        except Exception as e:
            print(f"Could not get real sequence: {e}")

    for i, seq in enumerate(test_sequences):
        print(f"\n--- Sequence {i+1} (length: {len(seq)}) ---")
        print(f"Seq: {seq[:50]}{'...' if len(seq) > 50 else ''}")
        
        try:
            # Get coverage without filtering
            coverage = index.get_sequence_coverage(seq, cutoff=0, k=23)
            
            if len(coverage) == 0:
                print("  Sequence too short for 23-mer analysis")
                continue
                
            # Coverage statistics
            non_zero_positions = sum(1 for tf in coverage if tf > 0)
            total_positions = len(coverage)
            max_tf = max(coverage) if coverage else 0
            avg_tf = sum(coverage) / len(coverage) if coverage else 0
            
            print(f"  Total positions for analysis: {total_positions}")
            print(f"  Positions with nonzero TF: {non_zero_positions} ({non_zero_positions/total_positions*100:.1f}%)")
            print(f"  Max TF: {max_tf}")
            print(f"  Average TF: {avg_tf:.2f}")
            
            # Test with different cutoffs
            cutoffs = [1, 5, 10] if max_tf > 0 else [0]
            for cutoff in cutoffs:
                filtered_coverage = index.get_sequence_coverage(seq, cutoff=cutoff, k=23)
                filtered_positions = sum(1 for tf in filtered_coverage if tf > 0)
                print(f"  With cutoff >= {cutoff}: {filtered_positions} positions")
            
            # Show first few k-mers with their TF
            print(f"  First 5 k-mers:")
            for j in range(min(5, len(coverage))):
                kmer = seq[j:j+23]
                tf = coverage[j]
                print(f"    Pos {j}: {kmer} -> TF = {tf}")
                
        except Exception as e:
            print(f"  Error analyzing coverage: {e}")

    # Test iterating over k-mers of a sequence
    print(f"\n=== Test iterating over sequence k-mers ===")
    
    # Instead of using test sequence, take a real read from index
    if index.n_reads > 0:
        # Get a random read from index
        rid = min(10, index.n_reads - 1)  # Take one of the first reads (not more than 10th)
        seq = index.get_read_by_rid(rid)
        
        # If read contains subreads, take first subread
        if '~' in seq:
            seq = seq.split('~')[0]
            
        # Limit length for demonstration
        if len(seq) > 200:
            seq = seq[:200]
            
        print(f"Analyzing read {rid}, length {len(seq)}")
        print(f"Read sequence: {seq[:50]}..." if len(seq) > 50 else f"Read sequence: {seq}")
        
        try:
            kmer_count = 0
            tf_sum = 0
            non_zero_kmers = 0
            
            for kmer, tf in index.iter_sequence_kmers(seq, k=23):
                kmer_count += 1
                tf_sum += tf
                if tf > 0:
                    non_zero_kmers += 1
                
                # Show first 5 k-mers with nonzero TF
                if tf > 0 and non_zero_kmers <= 5:
                    print(f"  K-mer {kmer_count}: {kmer} -> TF = {tf}")
                    
            if kmer_count > 0:
                coverage_percent = (non_zero_kmers / kmer_count) * 100 if kmer_count > 0 else 0
                print(f"Total k-mers: {kmer_count}")
                print(f"K-mers with nonzero TF: {non_zero_kmers} ({coverage_percent:.2f}%)")
                print(f"Average TF: {tf_sum/kmer_count:.2f}")
                
        except Exception as e:
            print(f"Error iterating over k-mers: {e}")
    else:
        print("Index contains no reads, skipping sequence k-mer iteration test")

    # Test print_sequence_coverage (if sequence is not too long)
    short_seq = test_sequences[0][:50] if len(test_sequences[0]) > 50 else test_sequences[0]
    if len(short_seq) >= 23:
        print(f"\n=== Test print_sequence_coverage for short sequence ===")
        print(f"Sequence: {short_seq}")
        
        try:
            print("Coverage:")
            coverage_result = index.print_sequence_coverage(short_seq, cutoff=0)
            print(f"Returned coverage array: {coverage_result}")
            
        except Exception as e:
            print(f"Error print_sequence_coverage: {e}")


# Test getting k-mer info by kid
if index is None:
    print("⚠️ Index not loaded, skipping k-mer info testing")
else:
    print("=== Test getting k-mer info by kid ===")

    if index.n_kmers == 0:
        print("⚠️ Index contains no k-mers")
    else:
        print(f"Total k-mers in index: {index.n_kmers:,}")
        
        # Test various kids
        test_kids = [0, 1, 10, 100, 1000, min(10000, index.n_kmers-1)]
        
        for kid in test_kids:
            if kid >= index.n_kmers:
                continue
                
            print(f"\n--- Kid {kid} ---")
            
            try:
                # Get k-mer by kid
                kmer = index.get_kmer_by_kid(kid)
                print(f"  Kmer: {kmer}")
                
                if kmer:
                    # Get TF for this k-mer
                    tf = index.get_tf_value(kmer)
                    print(f"  TF: {tf}")
                    
                    # Get hash value
                    hash_val = index.get_hash_value(kmer)
                    print(f"  Hash value: {hash_val}")
                    
                    # Get strand
                    strand = index.get_strand(kmer)
                    print(f"  Strand: {strand} ({'FORWARD' if strand == 1 else 'REVERSE' if strand == 2 else 'NOT_FOUND'})")
                    
                    # Get full info via get_kmer_info_by_kid
                    kmer_info, rkmer_info, tf_info = index.get_kmer_info_by_kid(kid)
                    print(f"  Info: kmer={kmer_info}, rkmer={rkmer_info}, tf={tf_info}")
                    
                    # Check consistency
                    if kmer == kmer_info and tf == tf_info:
                        print(f"  ✓ Data consistency confirmed")
                    else:
                        print(f"  ✗ Data mismatch!")
                        
                    # Reverse check: get kid by k-mer
                    kid_back = index.get_kid_by_kmer(kmer)
                    if kid_back == kid:
                        print(f"  ✓ Reverse conversion successful")
                    else:
                        print(f"  ✗ Reverse conversion failed: {kid} -> {kid_back}")
                        
                else:
                    print(f"  ✗ K-mer not found")
                    
            except Exception as e:
                print(f"  Error: {e}")

    # Test get_rid2poses function
    if index.aindex_loaded and 'test_kmers_for_positions' in globals() and test_kmers_for_positions:
        print(f"\n=== Test getting mapping RID -> positions ===")
        
        test_kmer, test_tf = test_kmers_for_positions[0]
        print(f"Testing k-mer: {test_kmer} (TF = {test_tf})")
        
        try:
            rid2poses = index.get_rid2poses(test_kmer)
            
            print(f"Reads with this k-mer: {len(rid2poses)}")
            
            # Show first few reads
            for i, (rid, poses) in enumerate(list(rid2poses.items())[:5]):
                print(f"  Read {rid}: positions in read = {poses}")
                
                # Check one of the positions
                if poses:
                    local_pos = poses[0]
                    try:
                        read = index.get_read_by_rid(rid)
                        if read and local_pos + 23 <= len(read):
                            subseq = read[local_pos:local_pos+23]
                            print(f"    Checking position {local_pos}: {subseq}")
                            if subseq == test_kmer:
                                print(f"    ✓ Match")
                            else:
                                print(f"    ✗ Mismatch")
                    except Exception as e:
                        print(f"    Error checking: {e}")
                        
        except Exception as e:
            print(f"Error get_rid2poses: {e}")

    # Final test summary
    print(f"\n" + "="*60)
    print(f"INDEX TESTING SUMMARY")
    print(f"="*60)
    prefix_display = f"{PATH_TO_AINDEX_FOLDER}/temp/reads.23"
    print(f"Prefix: {prefix_display}")
    print(f"Index loaded: {'✓' if index._loaded else '✗'}")
    print(f"AIndex loaded: {'✓' if index.aindex_loaded else '✗'}")
    print(f"Number of k-mers: {index.n_kmers:,}")
    print(f"Number of reads: {index.n_reads:,}")
    print(f"Reads file size: {index.reads_size:,} bytes")

    # Additional index info
    try:
        info = index.get_index_info()
        print(f"\nAdditional info:")
        print(info)
    except:
        pass

    print(f"\n✓ Testing completed successfully!")

# STRESS TEST: 1 MILLION QUERIES FOR 23-MERS
if index is None:
    print("⚠️ Index not loaded, skipping stress test")
else:
    print(f"\n" + "="*80)
    print(f"STRESS TEST: 1 MILLION QUERIES FOR 23-MERS")
    print(f"="*80)
    
    # Check for new functions for 23-mers (similar to 13-mers)
    has_23mer_functions = hasattr(index, 'get_tf_both_directions_23mer') and hasattr(index, 'get_total_tf_value_23mer')
    
    if not has_23mer_functions:
        print("⚠️ New functions for 23-mers not found, using standard functions")
    else:
        print("✓ New functions for 23-mers detected")
    
    # Generate 1 million random 23-mers
    NUM_TESTS = 1_000_000
    print(f"Generating {NUM_TESTS:,} random 23-mers...")
    
    start_gen = time.time()
    random_23mers = []
    nucleotides = ['A', 'T', 'G', 'C']
    
    for i in range(NUM_TESTS):
        kmer = ''.join(np.random.choice(nucleotides, 23))
        random_23mers.append(kmer)
        
        if (i + 1) % 100000 == 0:
            print(f"  Generated: {i+1:,} / {NUM_TESTS:,}")
    
    gen_time = time.time() - start_gen
    print(f"✓ Generation completed in {gen_time:.2f} seconds")
    
    # Test 1: Standard get_tf_value function (single queries)
    print(f"\n--- Test 1: Single queries get_tf_value ---")
    start_time = time.time()
    
    single_results = []
    for i, kmer in enumerate(random_23mers):
        try:
            tf = index.get_tf_value(kmer)
            single_results.append(tf)
        except:
            single_results.append(0)
            
        if (i + 1) % 100000 == 0:
            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed
            print(f"  Processed: {i+1:,} / {NUM_TESTS:,} ({rate:,.0f} queries/sec)")
    
    single_time = time.time() - start_time
    single_rate = NUM_TESTS / single_time
    non_zero_single = sum(1 for tf in single_results if tf > 0)
    
    print(f"✓ Single queries: {single_time:.2f} sec, {single_rate:,.0f} queries/sec")
    print(f"  Found k-mers with TF > 0: {non_zero_single:,} ({non_zero_single/NUM_TESTS*100:.2f}%)")
    
    # Test 2: Batch get_tf_values function
    print(f"\n--- Test 2: Batch queries get_tf_values ---")
    start_time = time.time()
    
    try:
        batch_results = index.get_tf_values(random_23mers)
        batch_time = time.time() - start_time
        batch_rate = NUM_TESTS / batch_time
        non_zero_batch = sum(1 for tf in batch_results if tf > 0)
        
        print(f"✓ Batch queries: {batch_time:.2f} sec, {batch_rate:,.0f} queries/sec")
        print(f"  Found k-mers with TF > 0: {non_zero_batch:,} ({non_zero_batch/NUM_TESTS*100:.2f}%)")
        print(f"  Speedup: {single_rate/batch_rate:.1f}x {'(batch faster)' if batch_rate > single_rate else '(single faster)'}")
        
        # Check consistency of results
        mismatches = sum(1 for i in range(len(single_results)) if single_results[i] != batch_results[i])
        print(f"  Consistency: {NUM_TESTS - mismatches:,} / {NUM_TESTS:,} matches ({(NUM_TESTS-mismatches)/NUM_TESTS*100:.3f}%)")
        
    except Exception as e:
        print(f"✗ Batch query error: {e}")
        batch_results = None
        batch_time = float('inf')
    
    # Test 3: New functions for 23-mers (if available)
    if has_23mer_functions:
        print(f"\n--- Test 3: New functions for 23-mers ---")
        
        # Test get_total_tf_value_23mer (single)
        print("Testing get_total_tf_value_23mer...")
        start_time = time.time()
        
        new_single_results = []
        for i, kmer in enumerate(random_23mers):
            try:
                tf = index.get_total_tf_value_23mer(kmer)
                new_single_results.append(tf)
            except:
                new_single_results.append(0)
                
            if (i + 1) % 100000 == 0:
                elapsed = time.time() - start_time
                rate = (i + 1) / elapsed
                print(f"  Processed: {i+1:,} / {NUM_TESTS:,} ({rate:,.0f} queries/sec)")
        
        new_single_time = time.time() - start_time
        new_single_rate = NUM_TESTS / new_single_time
        non_zero_new_single = sum(1 for tf in new_single_results if tf > 0)
        
        print(f"✓ get_total_tf_value_23mer: {new_single_time:.2f} sec, {new_single_rate:,.0f} queries/sec")
        print(f"  Found k-mers with TF > 0: {non_zero_new_single:,} ({non_zero_new_single/NUM_TESTS*100:.2f}%)")
        
        # Test get_total_tf_values_23mer (batch)
        if hasattr(index, 'get_total_tf_values_23mer'):
            print("Testing get_total_tf_values_23mer...")
            start_time = time.time()
            
            try:
                new_batch_results = index.get_total_tf_values_23mer(random_23mers)
                new_batch_time = time.time() - start_time
                new_batch_rate = NUM_TESTS / new_batch_time
                non_zero_new_batch = sum(1 for tf in new_batch_results if tf > 0)
                
                print(f"✓ get_total_tf_values_23mer: {new_batch_time:.2f} sec, {new_batch_rate:,.0f} queries/sec")
                print(f"  Found k-mers with TF > 0: {non_zero_new_batch:,} ({non_zero_new_batch/NUM_TESTS*100:.2f}%)")
                print(f"  Speedup: {new_single_rate/new_batch_rate:.1f}x {'(batch faster)' if new_batch_rate > new_single_rate else '(single faster)'}")
                
                # Compare with original functions
                if batch_results is not None:
                    diff_count = sum(1 for i in range(len(batch_results)) if batch_results[i] != new_batch_results[i])
                    print(f"  Differences with get_tf_values: {diff_count:,} ({diff_count/NUM_TESTS*100:.3f}%)")
                
            except Exception as e:
                print(f"✗ Error get_total_tf_values_23mer: {e}")
    
    # Performance analysis and statistics
    print(f"\n--- Statistics and performance analysis ---")
    
    # Stats for found k-mers
    if single_results:
        tf_values = [tf for tf in single_results if tf > 0]
        if tf_values:
            print(f"TF statistics for found k-mers:")
            print(f"  Min: {min(tf_values)}")
            print(f"  Max: {max(tf_values)}")
            print(f"  Mean: {np.mean(tf_values):.2f}")
            print(f"  Median: {np.median(tf_values):.2f}")
            print(f"  Std: {np.std(tf_values):.2f}")
        
    # Performance comparison
    print(f"\nPerformance ranking:")
    performance_data = [
        ("get_tf_value (single)", single_rate, single_time),
    ]
    
    if 'batch_rate' in locals() and batch_time != float('inf'):
        performance_data.append(("get_tf_values (batch)", batch_rate, batch_time))
    
    if has_23mer_functions and 'new_single_rate' in locals():
        performance_data.append(("get_total_tf_value_23mer", new_single_rate, new_single_time))
        
    if has_23mer_functions and 'new_batch_rate' in locals():
        performance_data.append(("get_total_tf_values_23mer", new_batch_rate, new_batch_time))
    
    # Sort by speed
    performance_data.sort(key=lambda x: x[1], reverse=True)
    
    print("Speed ranking:")
    for i, (name, rate, time_taken) in enumerate(performance_data, 1):
        print(f"  {i}. {name}: {rate:,.0f} queries/sec ({time_taken:.2f} sec)")
    
    # Memory usage (approximate)
    import psutil
    import os
    
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    print(f"\nMemory usage:")
    print(f"  RSS: {memory_info.rss / 1024 / 1024:.1f} MB")
    print(f"  VMS: {memory_info.vms / 1024 / 1024:.1f} MB")
    
    # Get index statistics (if available)
    if hasattr(index, 'get_23mer_statistics'):
        try:
            stats = index.get_23mer_statistics()
            print(f"\n23-mer index statistics:")
            for key, value in stats.items():
                if isinstance(value, (int, float)):
                    if isinstance(value, int) and value > 1000:
                        print(f"  {key}: {value:,}")
                    else:
                        print(f"  {key}: {value}")
                else:
                    print(f"  {key}: {value}")
        except Exception as e:
            print(f"Error getting statistics: {e}")
    
    print(f"\n" + "="*80)
    print(f"STRESS TEST COMPLETED")
    print(f"="*80)
    print(f"✓ Successfully tested {NUM_TESTS:,} queries for 23-mers")
    print(f"✓ Best performance: {max(p[1] for p in performance_data):,.0f} queries/sec")
    print(f"✓ Unique k-mers found: {non_zero_single:,} ({non_zero_single/NUM_TESTS*100:.2f}%)")

# Stress test for sequence coverage analysis
if index is not None:
    print("\n" + "="*80)
    print("STRESS TEST: COVERAGE ANALYSIS FOR 10,000 SEQUENCES")
    print("="*80)
    
    # Collect real sequences from index
    test_sequences = []
    print("Collecting test sequences from reads...")
    
    # Take first 10K reads (or as many as available)
    max_sequences = min(10000, index.n_reads)
    
    for rid in range(max_sequences):
        try:
            real_read = index.get_read_by_rid(rid)
            if real_read and len(real_read) >= 50:  # Minimum length for analysis
                # If read contains subreads, take first subread
                if '~' in real_read:
                    seq = real_read.split('~')[0]
                else:
                    seq = real_read
                
                # Limit length for stable runtime
                if len(seq) > 200:
                    seq = seq[:200]
                
                if len(seq) >= 23:  # Minimum length for 23-mers
                    test_sequences.append(seq)
                    
            # Show progress every 1000 sequences
            if (rid + 1) % 1000 == 0:
                print(f"  Collected: {len(test_sequences)} sequences from {rid + 1} reads")
                
        except Exception as e:
            continue  # Skip problematic reads
    
    print(f"✓ Collected {len(test_sequences)} sequences for testing")
    
    if len(test_sequences) == 0:
        print("⚠️ Could not collect sequences for testing")
    else:
        # Sequence statistics
        seq_lengths = [len(seq) for seq in test_sequences]
        avg_length = sum(seq_lengths) / len(seq_lengths)
        min_length = min(seq_lengths)
        max_length = max(seq_lengths)
        
        print(f"Sequence statistics:")
        print(f"  Count: {len(test_sequences):,}")
        print(f"  Average length: {avg_length:.1f} nt")
        print(f"  Length range: {min_length} - {max_length} nt")
        
        # Estimated number of k-mers
        total_kmers = sum(max(0, len(seq) - 22) for seq in test_sequences)
        print(f"  Total k-mers for analysis: {total_kmers:,}")
        
        print(f"\n--- Stress test coverage analysis ---")
        
        # Measure memory usage before test
        import psutil
        import os
        process = psutil.Process(os.getpid())
        memory_before = process.memory_info()
        
        print(f"Memory before test: RSS = {memory_before.rss / 1024 / 1024:.1f} MB")
        
        # Run test
        start_time = time.time()
        processed_sequences = 0
        total_positions_analyzed = 0
        total_nonzero_positions = 0
        total_tf_sum = 0
        max_tf_found = 0
        
        print("Analyzing sequence coverage...")
        
        for i, seq in enumerate(test_sequences):
            try:
                # Analyze coverage
                coverage = index.get_sequence_coverage(seq, cutoff=0, k=23)
                
                # Collect statistics
                processed_sequences += 1
                total_positions_analyzed += len(coverage)
                nonzero_positions = sum(1 for tf in coverage if tf > 0)
                total_nonzero_positions += nonzero_positions
                tf_sum = sum(coverage)
                total_tf_sum += tf_sum
                max_tf_in_seq = max(coverage) if coverage else 0
                max_tf_found = max(max_tf_found, max_tf_in_seq)
                
                # Show progress every 1000 sequences
                if processed_sequences % 1000 == 0:
                    elapsed = time.time() - start_time
                    sequences_per_sec = processed_sequences / elapsed if elapsed > 0 else 0
                    positions_per_sec = total_positions_analyzed / elapsed if elapsed > 0 else 0
                    
                    print(f"  Processed: {processed_sequences:,} / {len(test_sequences):,} " +
                          f"({sequences_per_sec:.1f} seq/sec, {positions_per_sec:.0f} positions/sec)")
                    
            except Exception as e:
                print(f"  Error analyzing sequence {i}: {e}")
                continue
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        
        # Measure memory after test
        memory_after = process.memory_info()
        memory_increase = (memory_after.rss - memory_before.rss) / 1024 / 1024
        
        print(f"\n--- Stress test results ---")
        print(f"✓ Sequences processed: {processed_sequences:,}")
        print(f"✓ Total runtime: {elapsed_time:.2f} seconds")
        print(f"✓ Performance:")
        print(f"  - Sequences per second: {processed_sequences / elapsed_time:.1f}")
        print(f"  - k-mer positions per second: {total_positions_analyzed / elapsed_time:.0f}")
        print(f"  - Time per sequence: {elapsed_time * 1000 / processed_sequences:.2f} ms")
        
        print(f"\n✓ Coverage statistics:")
        print(f"  - Total positions analyzed: {total_positions_analyzed:,}")
        print(f"  - Positions with nonzero TF: {total_nonzero_positions:,} ({total_nonzero_positions/total_positions_analyzed*100:.1f}%)")
        print(f"  - Average TF per position: {total_tf_sum/total_positions_analyzed:.2f}")
        print(f"  - Max TF: {max_tf_found}")
        print(f"  - Total TF sum: {total_tf_sum:,}")
        
        print(f"\n✓ Memory usage:")
        print(f"  - Memory after test: RSS = {memory_after.rss / 1024 / 1024:.1f} MB")
        print(f"  - Memory increase: {memory_increase:+.1f} MB")
        
        # Additional test: several sequences with different cutoffs
        print(f"\n--- Test effect of cutoff on performance ---")
        test_sample = test_sequences[:100]  # Take 100 sequences
        cutoff_values = [0, 1, 5, 10, 20]
        
        for cutoff in cutoff_values:
            start_cutoff_time = time.time()
            total_positions_cutoff = 0
            
            for seq in test_sample:
                try:
                    coverage = index.get_sequence_coverage(seq, cutoff=cutoff, k=23)
                    total_positions_cutoff += sum(1 for tf in coverage if tf > 0)
                except:
                    continue
            
            cutoff_time = time.time() - start_cutoff_time
            print(f"  Cutoff >= {cutoff}: {cutoff_time:.3f} sec, {total_positions_cutoff} positions with TF")
        
        print(f"\n" + "="*80)
        print(f"COVERAGE STRESS TEST COMPLETED")
        print(f"="*80)
        print(f"✓ Successfully analyzed {processed_sequences:,} sequences")
        print(f"✓ Performance: {processed_sequences / elapsed_time:.1f} sequences/sec")
        print(f"✓ Total time: {elapsed_time:.2f} seconds")
else:
    print("\n⚠️ Index not loaded, skipping coverage stress test")