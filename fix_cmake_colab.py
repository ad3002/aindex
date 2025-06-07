#!/usr/bin/env python3
"""
Quick fix for cmake conflict in Google Colab
Run this before installing aindex
"""

import subprocess
import sys

def run_cmd(cmd):
    print(f"$ {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {e}")
        if e.stderr:
            print("STDERR:", e.stderr)
        return False

print("=== FIXING CMAKE CONFLICT IN GOOGLE COLAB ===\n")

# Check if Python cmake is installed
try:
    import cmake
    print("❌ Found Python cmake module - this will cause conflicts!")
    
    print("\n1. Removing Python cmake module...")
    if not run_cmd("pip uninstall -y cmake"):
        print("Failed to remove Python cmake")
        
    print("\n2. Installing system cmake...")
    if not run_cmd("apt-get update"):
        print("Failed to update package lists")
    
    if not run_cmd("apt-get install -y cmake"):
        print("Failed to install system cmake")
        
    print("\n3. Verifying cmake installation...")
    run_cmd("which cmake")
    run_cmd("cmake --version")
    
    print("\n✅ CMAKE CONFLICT RESOLVED!")
    print("Now you can install aindex with: pip install .")
    
except ImportError:
    print("✅ No Python cmake module found - no conflict to resolve")
    print("Ensuring system cmake is installed...")
    run_cmd("apt-get update && apt-get install -y cmake")
    run_cmd("cmake --version")

print("\n=== QUICK FIX COMPLETE ===")
