#!/usr/bin/env python3
"""
Test script to verify input/output validation in aindex CLI
"""

import tempfile
import os
import sys
from pathlib import Path

# Add aindex to path for testing
sys.path.insert(0, str(Path(__file__).parent))

from aindex.cli import validate_input_output_files, validate_output_file_overwrite

def test_validation():
    """Test the validation functions"""
    
    print("=== Testing Input/Output Validation ===")
    
    # Test 1: Same input and output files
    print("\nTest 1: Same input and output files")
    result = validate_input_output_files("test.fastq", "test.fastq", "test-command")
    print(f"Result: {result} (should be False)")
    
    # Test 2: Different input and output files
    print("\nTest 2: Different input and output files")
    result = validate_input_output_files("input.fastq", "output.txt", "test-command")
    print(f"Result: {result} (should be True)")
    
    # Test 3: Similar names (warning case)
    print("\nTest 3: Similar names")
    result = validate_input_output_files("data.fastq", "data.fastq.counts", "test-command")
    print(f"Result: {result} (should be True with warning)")
    
    # Test 4: Empty files
    print("\nTest 4: Empty files")
    result = validate_input_output_files("", "output.txt", "test-command")
    print(f"Result: {result} (should be True)")
    
    # Test 5: Test output file overwrite with existing file
    print("\nTest 5: Output file overwrite validation")
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
        tmp.write("test content")
        tmp_path = tmp.name
    
    try:
        result = validate_output_file_overwrite(tmp_path, "test-command")
        print(f"Result: {result} (should be True with warning)")
    finally:
        os.unlink(tmp_path)
    
    # Test 6: Non-existent output file
    print("\nTest 6: Non-existent output file")
    result = validate_output_file_overwrite("/tmp/non_existent_file.txt", "test-command")
    print(f"Result: {result} (should be True)")
    
    print("\n=== Tests completed ===")

if __name__ == "__main__":
    test_validation()
