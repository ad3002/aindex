#!/usr/bin/env python3
"""
Installation script for Google Colab environment
Usage in Colab:
    !wget https://raw.githubusercontent.com/ad3002/aindex/main/install_colab.py
    !python install_colab.py
"""

import subprocess
import sys
import os

def run_command(cmd, description="", capture_output=True, check=True):
    """Run a command and handle errors gracefully"""
    print(f"\n{'='*60}")
    print(f"Running: {description or cmd}")
    print(f"Command: {cmd}")
    print(f"{'='*60}")
    
    try:
        if capture_output:
            result = subprocess.run(cmd, shell=True, check=check, 
                                  capture_output=True, text=True)
            if result.stdout:
                print("STDOUT:")
                print(result.stdout)
            if result.stderr:
                print("STDERR:")
                print(result.stderr)
            return result.returncode == 0
        else:
            result = subprocess.run(cmd, shell=True, check=check)
            return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code: {e.returncode}")
        if hasattr(e, 'stdout') and e.stdout:
            print("STDOUT:")
            print(e.stdout)
        if hasattr(e, 'stderr') and e.stderr:
            print("STDERR:")
            print(e.stderr)
        return False

def debug_build_failure():
    """Try to diagnose build failure step by step"""
    print("\n" + "="*60)
    print("DEBUGGING BUILD FAILURE")
    print("="*60)
    
    # Check Python development headers
    print("\n1. Checking Python development headers...")
    run_command("python3 -c 'import sysconfig; print(sysconfig.get_path(\"include\"))'", 
                "Python include path")
    
    # Check if we can compile a simple C++ program
    print("\n2. Testing C++ compilation...")
    test_cpp = """
#include <iostream>
int main() {
    std::cout << "Hello World" << std::endl;
    return 0;
}
"""
    with open("test.cpp", "w") as f:
        f.write(test_cpp)
    
    if run_command("g++ -o test test.cpp", "Simple C++ compilation test"):
        run_command("./test", "Running test program")
    run_command("rm -f test test.cpp", "Cleanup test files")
    
    # Check Python extension compilation
    print("\n3. Testing Python extension compilation...")
    test_pyx = """
#include <Python.h>

static PyObject* hello(PyObject* self, PyObject* args) {
    return PyUnicode_FromString("Hello from C++!");
}

static PyMethodDef methods[] = {
    {"hello", hello, METH_NOARGS, "Say hello"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "test_ext",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC PyInit_test_ext(void) {
    return PyModule_Create(&module);
}
"""
    with open("test_ext.cpp", "w") as f:
        f.write(test_pyx)
    
    # Get Python include directory
    python_include = subprocess.run("python3 -c 'import sysconfig; print(sysconfig.get_path(\"include\"))'", 
                                   shell=True, capture_output=True, text=True).stdout.strip()
    
    compile_cmd = f"g++ -shared -fPIC -I{python_include} test_ext.cpp -o test_ext.so"
    if run_command(compile_cmd, "Python extension compilation test"):
        print("✓ Python extension compilation works")
    else:
        print("✗ Python extension compilation failed")
        # Try to install python3-dev
        print("Installing python3-dev...")
        run_command("apt-get install -y python3-dev", "Installing Python development headers")
        run_command(compile_cmd, "Retry Python extension compilation")
    
    run_command("rm -f test_ext.cpp test_ext.so", "Cleanup extension test files")
    
    # Try manual build
    print("\n4. Attempting manual build...")
    if os.path.exists("Makefile"):
        print("Found Makefile, trying manual build...")
        run_command("make clean", "Clean previous build", capture_output=False, check=False)
        if not run_command("make external", "Build external dependencies", capture_output=False, check=False):
            print("External dependencies build failed. Checking git access...")
            run_command("git config --global http.timeout 300", "Increase git timeout")
            run_command("git config --global http.postBuffer 1048576000", "Increase git buffer")
            run_command("make external", "Retry building external dependencies", capture_output=False, check=False)
        
        run_command("make all", "Build all components", capture_output=False, check=False)
    
    print("\n5. Checking build artifacts...")
    run_command("find . -name '*.so' -o -name '*.exe'", "Find built binaries")

def install_colab():
    """Install aindex in Google Colab"""
    
    print("=== Installing aindex in Google Colab ===\n")
    
    # Install system dependencies
    print("1. Installing system dependencies...")
    if not run_command("apt-get update", "Updating package lists"):
        return False
    
    if not run_command("apt-get install -y build-essential cmake git python3-dev", 
                      "Installing build tools"):
        return False
    
    # Clone the repository if not exists
    if not os.path.exists("aindex"):
        print("\n2. Cloning aindex repository...")
        if not run_command("git clone https://github.com/ad3002/aindex.git", 
                          "Cloning repository"):
            return False
        os.chdir("aindex")
    else:
        print("\n2. Using existing aindex directory...")
        os.chdir("aindex")
    
    # Install Python dependencies
    print("\n3. Installing Python dependencies...")
    if not run_command("pip install -r requirements.txt", 
                      "Installing Python requirements"):
        return False
    
    # Try to install the package
    print("\n4. Installing the package...")
    if run_command("pip install -e . --verbose", 
                   "Installing aindex package with verbose output", 
                   capture_output=False, check=False):
        print("\n=== Installation completed successfully! ===")
        return True
    else:
        print("\n=== Installation failed. Running diagnostics... ===")
        debug_build_failure()
        return False

if __name__ == "__main__":
    # Check if we're in Colab
    try:
        import google.colab
        print("Google Colab environment detected.")
    except ImportError:
        print("Warning: This script is designed for Google Colab.")
    
    success = install_colab()
    sys.exit(0 if success else 1)
