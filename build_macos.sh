#!/bin/bash
set -ex

# This script is called by cibuildwheel.
# The first argument ($1) is the target architecture, e.g., "x86_64" or "arm64".

# Install build dependencies for the correct target architecture.
# The ARCHFLAGS environment variable (set in pyproject.toml) ensures
# pip installs the correct wheel.
pip install pybind11 numpy

# Clean any previous build artifacts.
make clean || true

# Recreate necessary directories.
mkdir -p bin obj aindex/core aindex/bin

# Select the correct make target based on the architecture.
if [ "$1" = "arm64" ]; then
    echo "--- Building ARM64 version for Apple Silicon ---"
    make arm64
else
    echo "--- Building x86_64 version ---"
    make all
fi

echo "--- Build complete, listing artifacts for verification ---"
ls -la aindex/core/ || echo "aindex/core not found"
ls -la bin/ || echo "bin not found"