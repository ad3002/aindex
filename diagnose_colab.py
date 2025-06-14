#!/usr/bin/env python3
"""
Quick diagnostic script for aindex build issues in Google Colab
Run this to diagnose why the build is failing
"""

import subprocess
import sys
import os

def run_cmd(cmd):
    """Run command and show output"""
    print(f"\n$ {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)
        return result.returncode == 0
    except Exception as e:
        print(f"Error running command: {e}")
        return False

print("=== AINDEX BUILD DIAGNOSTICS ===\n")

# Check system info
print("1. System Information:")
run_cmd("uname -a")
run_cmd("cat /etc/os-release | head -5")

# Check build tools
print("\n2. Build Tools:")
run_cmd("which gcc")
run_cmd("gcc --version | head -1")
run_cmd("which g++")
run_cmd("g++ --version | head -1")
run_cmd("which make")
run_cmd("make --version | head -1")

# Check cmake - this is often problematic in Colab
print("\n2.1. CMake Analysis:")
run_cmd("which cmake")
run_cmd("ls -la /usr/bin/cmake || echo 'System cmake not found'")
run_cmd("ls -la /usr/local/bin/cmake || echo 'Local cmake not found'")
if run_cmd("python3 -c 'import cmake' 2>/dev/null"):
    print("⚠️  WARNING: Python cmake module is installed!")
    print("This can conflict with system cmake. Consider:")
    print("pip uninstall cmake && apt-get install cmake")
else:
    print("✓ No Python cmake module found")

run_cmd("cmake --version | head -1")
run_cmd("which git")
run_cmd("git --version")

# Check Python
print("\n3. Python Information:")
run_cmd("which python3")
run_cmd("python3 --version")
run_cmd("python3 -c 'import sys; print(sys.executable)'")
run_cmd("python3 -c 'import sysconfig; print(sysconfig.get_path(\"include\"))'")

# Check Python development headers
print("\n4. Python Development Headers:")
run_cmd("find /usr/include -name 'Python.h' 2>/dev/null | head -5")
run_cmd("pkg-config --exists python3 && echo 'python3.pc found' || echo 'python3.pc NOT found'")

# Test simple compilation
print("\n5. Compilation Test:")
test_cpp = '''#include <iostream>
#include <Python.h>
int main() { std::cout << "C++ and Python headers OK" << std::endl; return 0; }'''

with open("/tmp/test.cpp", "w") as f:
    f.write(test_cpp)

python_include = subprocess.run("python3 -c 'import sysconfig; print(sysconfig.get_path(\"include\"))'", 
                               shell=True, capture_output=True, text=True).stdout.strip()

if run_cmd(f"g++ -I{python_include} /tmp/test.cpp -o /tmp/test"):
    run_cmd("/tmp/test")
    print("✓ Basic compilation works")
else:
    print("✗ Basic compilation failed")
    print("Installing python3-dev...")
    run_cmd("apt-get update && apt-get install -y python3-dev")
    print("Retrying compilation...")
    run_cmd(f"g++ -I{python_include} /tmp/test.cpp -o /tmp/test")

# Check if in aindex directory and try manual build
print("\n6. Manual Build Test:")
if os.path.exists("Makefile"):
    print("Found Makefile in current directory")
    run_cmd("make clean")
    
    # Check for cmake conflicts first
    print("\n6.1. Resolving cmake conflicts...")
    if run_cmd("python3 -c 'import cmake' 2>/dev/null"):
        print("Found Python cmake module - this may cause conflicts!")
        print("Removing Python cmake and installing system cmake...")
        run_cmd("pip uninstall -y cmake")
        run_cmd("apt-get update && apt-get install -y cmake")
    
    # Try to build external dependencies first  
    print("\n6.2. Building external dependencies...")
    if not run_cmd("timeout 300 make external"):
        print("External build failed or timed out. Trying manual approach...")
        run_cmd("mkdir -p external bin aindex/core")
        
        if not os.path.exists("external/emphf"):
            run_cmd("cd external && git clone https://github.com/ad3002/emphf.git")
        
        print("Building emphf manually...")
        run_cmd("cd external/emphf && ls -la")
        if run_cmd("cd external/emphf && /usr/bin/cmake . && make"):
            run_cmd("cp external/emphf/compute_mphf_seq bin/")
            print("✓ emphf build successful")
        else:
            print("✗ emphf build failed")
    
    print("\n6.3. Building main components...")
    run_cmd("make all")
    
    print("\n6.4. Checking build artifacts...")
    run_cmd("find . -name '*.so' -o -name '*.exe' | head -10")
else:
    print("No Makefile found in current directory")
    run_cmd("ls -la")

print("\n7. Pip install with verbose output:")
print("Now trying: pip install -v .")
run_cmd("pip install -v .")

print("\n=== DIAGNOSTICS COMPLETE ===")
print("If you see any obvious failures above, fix them first.")
print("Otherwise, the verbose pip output should show the exact error.")
