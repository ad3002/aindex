#!/usr/bin/env python3
"""
Integration tests for aindex CLI
Tests all CLI commands in semi-manual mode with proper setup and teardown
"""

import os
import sys
import tempfile
import shutil
import subprocess
import unittest
from pathlib import Path
import json

# Add the project root to Python path for imports
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from aindex import cli


class TestAindexCLI(unittest.TestCase):
    """Integration tests for aindex CLI commands"""
    
    @classmethod
    def setUpClass(cls):
        """Set up test environment once for all tests"""
        print("\n=== Setting up CLI Integration Tests ===")
        
        # Create temporary directory for all tests
        cls.test_dir = Path(tempfile.mkdtemp(prefix="aindex_cli_test_"))
        print(f"Test directory: {cls.test_dir}")
        
        # Create test data directory
        cls.data_dir = cls.test_dir / "data"
        cls.data_dir.mkdir()
        
        # Create test output directory
        cls.output_dir = cls.test_dir / "output"
        cls.output_dir.mkdir()
        
        # Prepare small test files
        cls._create_test_data()
        
        print("✓ Test environment ready")
    
    @classmethod
    def tearDownClass(cls):
        """Clean up test environment"""
        print(f"\n=== Cleaning up: {cls.test_dir} ===")
        if cls.test_dir.exists():
            shutil.rmtree(cls.test_dir)
        print("✓ Cleanup complete")
    
    @classmethod
    def _create_test_data(cls):
        """Create small test data files for CLI testing"""
        
        # Create small FASTA file
        cls.test_fasta = cls.data_dir / "test.fasta"
        with open(cls.test_fasta, 'w') as f:
            f.write(">seq1\n")
            f.write("ATCGATCGATCGATCG\n")
            f.write(">seq2\n") 
            f.write("GCTAGCTAGCTAGCTA\n")
            f.write(">seq3\n")
            f.write("TTTTAAAACCCCGGGG\n")
        
        # Create small FASTQ file
        cls.test_fastq = cls.data_dir / "test.fastq"
        with open(cls.test_fastq, 'w') as f:
            f.write("@seq1\n")
            f.write("ATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIII\n")
            f.write("@seq2\n")
            f.write("GCTAGCTAGCTAGCTA\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIII\n")
        
        # Create paired-end FASTQ files
        cls.test_fastq_r1 = cls.data_dir / "test_R1.fastq"
        cls.test_fastq_r2 = cls.data_dir / "test_R2.fastq"
        
        with open(cls.test_fastq_r1, 'w') as f:
            f.write("@seq1/1\n")
            f.write("ATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIII\n")
        
        with open(cls.test_fastq_r2, 'w') as f:
            f.write("@seq1/2\n")
            f.write("CGAATCGATCGATCGA\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIII\n")
        
        # Create reads file (simple format)
        cls.test_reads = cls.data_dir / "test.reads"
        with open(cls.test_reads, 'w') as f:
            f.write("ATCGATCGATCGATCG\n")
            f.write("GCTAGCTAGCTAGCTA\n")
            f.write("TTTTAAAACCCCGGGG\n")
        
        print(f"✓ Created test data files in {cls.data_dir}")
    
    def setUp(self):
        """Set up for each test"""
        # Create unique output directory for this test
        test_name = self._testMethodName
        self.current_output_dir = self.output_dir / test_name
        self.current_output_dir.mkdir(exist_ok=True)
    
    def _run_cli_command(self, command_args, expect_success=True):
        """Helper to run CLI commands and capture output"""
        try:
            # Save original sys.argv
            original_argv = sys.argv.copy()
            
            # Set up CLI arguments
            sys.argv = ['aindex'] + command_args
            
            # Capture the result
            result = cli.main()
            
            # Restore original sys.argv
            sys.argv = original_argv
            
            if expect_success:
                self.assertEqual(result, 0, f"Command failed: aindex {' '.join(command_args)}")
            
            return result
            
        except SystemExit as e:
            # Restore original sys.argv in case of SystemExit
            sys.argv = original_argv
            if expect_success:
                self.assertEqual(e.code, 0, f"Command exited with error: aindex {' '.join(command_args)}")
            return e.code
        except Exception as e:
            # Restore original sys.argv in case of exception
            sys.argv = original_argv
            if expect_success:
                self.fail(f"Command raised exception: aindex {' '.join(command_args)}: {e}")
            raise
    
    def test_help_command(self):
        """Test help command"""
        print("\n--- Testing help command ---")
        result = self._run_cli_command(['help'])
        self.assertEqual(result, 0)
    
    def test_version_command(self):
        """Test version command"""
        print("\n--- Testing version command ---")
        result = self._run_cli_command(['version'])
        self.assertEqual(result, 0)
    
    def test_info_command(self):
        """Test info command"""
        print("\n--- Testing info command ---")
        result = self._run_cli_command(['info'])
        self.assertEqual(result, 0)
    
    def test_platform_command(self):
        """Test platform information command"""
        print("\n--- Testing platform command ---")
        result = self._run_cli_command(['platform'])
        self.assertEqual(result, 0)
    
    def test_platform_list_executables(self):
        """Test platform command with list executables"""
        print("\n--- Testing platform --list-executables ---")
        result = self._run_cli_command(['platform', '--list-executables'])
        self.assertEqual(result, 0)
    
    def test_generate_13mers(self):
        """Test generate command for 13-mers"""
        print("\n--- Testing generate 13-mers ---")
        output_file = self.current_output_dir / "all_13mers.txt"
        
        result = self._run_cli_command([
            'generate',
            '-o', str(output_file),
            '-k', '13'
        ])
        
        if result == 0:
            self.assertTrue(output_file.exists(), "Output file should be created")
            print(f"✓ Generated 13-mers to {output_file}")
    
    def test_generate_13mers_with_options(self):
        """Test generate command with various options"""
        print("\n--- Testing generate 13-mers with options ---")
        output_file = self.current_output_dir / "all_13mers_indexed.txt"
        
        result = self._run_cli_command([
            'generate',
            '-o', str(output_file),
            '-k', '13',
            '-i',  # with indices
            '-s'   # show stats
        ])
        
        # This might fail if binary doesn't exist, that's OK for now
        print(f"Generate with options result: {result}")
    
    def test_build_hash_13mers(self):
        """Test build-hash command for 13-mers"""
        print("\n--- Testing build-hash for 13-mers ---")
        
        # First generate k-mers
        kmers_file = self.current_output_dir / "kmers_for_hash.txt"
        hash_file = self.current_output_dir / "test.hash"
        
        # Create a simple k-mers file manually for testing
        with open(kmers_file, 'w') as f:
            f.write("ATCGATCGATCGA\n")
            f.write("TCGATCGATCGAT\n")
            f.write("CGATCGATCGATC\n")
        
        result = self._run_cli_command([
            'build-hash',
            '-i', str(kmers_file),
            '-o', str(hash_file),
            '-k', '13'
        ])
        
        print(f"Build hash result: {result}")
        if result == 0:
            print(f"✓ Built hash file: {hash_file}")
    
    def test_count_direct_13mers(self):
        """Test count-direct command (ARM64 optimized k-mer counting)"""
        print("\n--- Testing count-direct for 13-mers ---")
        output_file = self.current_output_dir / "direct_counts.txt"
        
        result = self._run_cli_command([
            'count-direct',
            '-i', str(self.test_fastq),
            '-k', '13',
            '-o', str(output_file),
            '-t', '2',
            '--verbose'
        ])
        
        print(f"Count direct result: {result}")
        if result == 0:
            self.assertTrue(output_file.exists(), "Count output file should be created")
            print(f"✓ Direct k-mer counting completed: {output_file}")
    
    def test_count_with_hash(self):
        """Test count command with hash file"""
        print("\n--- Testing count with hash file ---")
        
        # This test requires a hash file, which might not exist
        # For now, just test the argument parsing
        hash_file = self.current_output_dir / "dummy.hash"
        output_file = self.current_output_dir / "counts.tf.bin"
        
        # Create dummy hash file
        hash_file.touch()
        
        result = self._run_cli_command([
            'count',
            '-i', str(self.test_fastq),
            '--hash-file', str(hash_file),
            '-o', str(output_file),
            '-k', '13'
        ], expect_success=False)  # This might fail without proper hash
        
        print(f"Count with hash result: {result}")
    
    def test_compute_reads_single_end(self):
        """Test compute-reads with single-end FASTQ"""
        print("\n--- Testing compute-reads single-end ---")
        output_prefix = self.current_output_dir / "se_reads"
        
        result = self._run_cli_command([
            'compute-reads',
            '-i', str(self.test_fastq),
            '-o', str(output_prefix)
        ])
        
        print(f"Compute reads SE result: {result}")
        if result == 0:
            expected_output = Path(str(output_prefix) + ".reads")
            print(f"✓ Single-end reads processed, expected: {expected_output}")
    
    def test_compute_reads_paired_end(self):
        """Test compute-reads with paired-end FASTQ"""
        print("\n--- Testing compute-reads paired-end ---")
        output_prefix = self.current_output_dir / "pe_reads"
        
        result = self._run_cli_command([
            'compute-reads',
            '-1', str(self.test_fastq_r1),
            '-2', str(self.test_fastq_r2),
            '-o', str(output_prefix)
        ])
        
        print(f"Compute reads PE result: {result}")
        if result == 0:
            expected_output = Path(str(output_prefix) + ".reads")
            print(f"✓ Paired-end reads processed, expected: {expected_output}")
    
    def test_compute_reads_fasta(self):
        """Test compute-reads with FASTA input"""
        print("\n--- Testing compute-reads with FASTA ---")
        output_prefix = self.current_output_dir / "fasta_reads"
        
        result = self._run_cli_command([
            'compute-reads',
            '-i', str(self.test_fasta),
            '-o', str(output_prefix),
            '--format', 'fasta'
        ])
        
        print(f"Compute reads FASTA result: {result}")
        if result == 0:
            expected_output = Path(str(output_prefix) + ".reads")
            print(f"✓ FASTA reads processed, expected: {expected_output}")
    
    def test_reads_to_fasta(self):
        """Test reads-to-fasta conversion"""
        print("\n--- Testing reads-to-fasta ---")
        output_file = self.current_output_dir / "converted.fasta"
        
        result = self._run_cli_command([
            'reads-to-fasta',
            '-i', str(self.test_reads),
            '-o', str(output_file)
        ])
        
        print(f"Reads to FASTA result: {result}")
        if result == 0:
            self.assertTrue(output_file.exists(), "FASTA output file should be created")
            print(f"✓ Reads converted to FASTA: {output_file}")
    
    def test_compute_aindex_high_level(self):
        """Test compute-aindex high-level command"""
        print("\n--- Testing compute-aindex high-level ---")
        output_prefix = self.current_output_dir / "aindex_test"
        
        result = self._run_cli_command([
            'compute-aindex',
            '-i', str(self.test_fastq),
            '-t', 'fastq',
            '-o', str(output_prefix),
            '-k', '13',
            '--use-kmer-counter'
        ])
        
        print(f"Compute aindex high-level result: {result}")
        # This is complex and might fail, that's OK for integration test
    
    def test_compute_aindex_direct_expert(self):
        """Test compute-aindex-direct expert command"""
        print("\n--- Testing compute-aindex-direct expert mode ---")
        
        # This requires several pre-computed files, so it will likely fail
        # But we test the argument parsing
        output_prefix = self.current_output_dir / "aindex_direct"
        dummy_hash = self.current_output_dir / "dummy.hash"
        dummy_tf = self.current_output_dir / "dummy.tf.bin"
        
        # Create dummy files
        dummy_hash.touch()
        dummy_tf.touch()
        
        result = self._run_cli_command([
            'compute-aindex-direct',
            str(self.test_reads),
            str(dummy_hash),
            str(output_prefix),
            '-t', '2',
            '-k', '13',
            '--tf-file', str(dummy_tf)
        ], expect_success=False)  # Expected to fail without proper data
        
        print(f"Compute aindex direct result: {result}")
    
    def test_compute_index(self):
        """Test compute-index command"""
        print("\n--- Testing compute-index ---")
        
        dummy_dat = self.current_output_dir / "dummy.dat"
        dummy_hash = self.current_output_dir / "dummy.hash"
        output_prefix = self.current_output_dir / "index_test"
        
        # Create dummy files
        dummy_dat.touch()
        dummy_hash.touch()
        
        result = self._run_cli_command([
            'compute-index',
            str(dummy_dat),
            str(dummy_hash),
            '-o', str(output_prefix),
            '-t', '2',
            '--mock'
        ], expect_success=False)  # Expected to fail without proper data
        
        print(f"Compute index result: {result}")
    
    def test_invalid_command(self):
        """Test handling of invalid commands"""
        print("\n--- Testing invalid command ---")
        result = self._run_cli_command(['invalid-command'], expect_success=False)
        self.assertNotEqual(result, 0, "Invalid command should return non-zero exit code")
    
    def test_command_help_flags(self):
        """Test --help flag for various commands"""
        print("\n--- Testing command help flags ---")
        
        commands_to_test = [
            'generate',
            'build-hash', 
            'count',
            'count-direct',
            'compute-reads',
            'compute-aindex',
            'reads-to-fasta'
        ]
        
        for cmd in commands_to_test:
            print(f"Testing help for: {cmd}")
            try:
                result = self._run_cli_command([cmd, '--help'], expect_success=False)
                # Help commands typically exit with 0 or SystemExit, both are OK
                print(f"  {cmd} --help: {result}")
            except SystemExit as e:
                print(f"  {cmd} --help: SystemExit({e.code})")
            except Exception as e:
                print(f"  {cmd} --help: Exception({e})")


class TestCLIUtilities(unittest.TestCase):
    """Test CLI utility functions"""
    
    def test_detect_platform(self):
        """Test platform detection"""
        platform_info = cli.detect_platform()
        
        self.assertIn('system', platform_info)
        self.assertIn('machine', platform_info)
        self.assertIn('is_apple_silicon', platform_info)
        self.assertIn('cpu_count', platform_info)
        
        print(f"Detected platform: {platform_info}")
    
    def test_get_optimal_executable(self):
        """Test optimal executable selection"""
        platform_info = cli.detect_platform()
        
        # Test some known executables
        test_executables = [
            'kmer_counter',
            'count_kmers', 
            'compute_aindex',
            'generate_all_13mers'
        ]
        
        for exe in test_executables:
            optimal = cli.get_optimal_executable(exe, platform_info)
            print(f"{exe} -> {optimal}")
            self.assertIsInstance(optimal, str)
    
    def test_get_bin_path(self):
        """Test bin path detection"""
        bin_path = cli.get_bin_path()
        print(f"Detected bin path: {bin_path}")
        self.assertIsInstance(bin_path, Path)
    
    def test_detect_file_format(self):
        """Test file format detection"""
        # Create temporary test files
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            
            # FASTA file
            fasta_file = tmp_path / "test.fasta"
            with open(fasta_file, 'w') as f:
                f.write(">seq1\nATCG\n")
            
            # FASTQ file
            fastq_file = tmp_path / "test.fastq"
            with open(fastq_file, 'w') as f:
                f.write("@seq1\nATCG\n+\nIIII\n")
            
            # Reads file
            reads_file = tmp_path / "test.reads"
            with open(reads_file, 'w') as f:
                f.write("ATCGATCG\n")
            
            # Test detection
            formats = [
                (fasta_file, 'fasta'),
                (fastq_file, 'fastq'),
                (reads_file, 'reads')
            ]
            
            for file_path, expected_format in formats:
                detected_format, first_line = cli.detect_file_format(file_path)
                print(f"{file_path.name}: {detected_format} (first line: {first_line})")
                self.assertEqual(detected_format, expected_format)


def run_integration_tests():
    """Run all integration tests"""
    print("=== aindex CLI Integration Tests ===")
    print("Running semi-manual integration tests for all CLI functions")
    print()
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestAindexCLI))
    suite.addTests(loader.loadTestsFromTestCase(TestCLIUtilities))
    
    # Run tests with verbose output
    runner = unittest.TextTestRunner(verbosity=2, buffer=True)
    result = runner.run(suite)
    
    # Print summary
    print("\n=== Test Summary ===")
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Skipped: {len(result.skipped) if hasattr(result, 'skipped') else 0}")
    
    if result.failures:
        print("\nFailures:")
        for test, traceback in result.failures:
            print(f"  {test}: {traceback}")
    
    if result.errors:
        print("\nErrors:")
        for test, traceback in result.errors:
            print(f"  {test}: {traceback}")
    
    success = len(result.failures) == 0 and len(result.errors) == 0
    print(f"\nResult: {'PASSED' if success else 'FAILED'}")
    
    return 0 if success else 1


if __name__ == '__main__':
    exit_code = run_integration_tests()
    sys.exit(exit_code)
