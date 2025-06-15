#!/bin/bash

# aindex-arm64: Wrapper script for ARM64-optimized aindex tools
# Automatically selects the best tools for Apple Silicon

# Check if we're on ARM64 Mac
if [[ $(uname -s) == "Darwin" && $(uname -m) == "arm64" ]]; then
    echo "🚀 ARM64-optimized aindex tools for Apple Silicon"
    echo "Platform: $(uname -s) $(uname -m)"
    echo "CPU cores: $(sysctl -n hw.ncpu)"
    echo
    
    # Set default thread count for M1/M2
    DEFAULT_THREADS=$(sysctl -n hw.ncpu)
    export AINDEX_THREADS=${AINDEX_THREADS:-$DEFAULT_THREADS}
    
    echo "Using $AINDEX_THREADS threads for optimal performance"
    echo
    
    # Run the CLI with ARM64 optimizations
    python3 "$(dirname "$0")/aindex/cli.py" "$@"
else
    echo "Note: This script is optimized for ARM64 Macs"
    echo "Platform: $(uname -s) $(uname -m)"
    echo "Using standard aindex tools"
    echo
    
    # Run the CLI with standard tools
    python3 "$(dirname "$0")/aindex/cli.py" "$@"
fi
