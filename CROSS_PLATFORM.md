# Cross-Platform Compatibility Guide

## Overview
This document describes the cross-platform compatibility approach for the aindex project to ensure it builds and runs correctly on both macOS M1 (ARM64) and Ubuntu (x86_64) systems.

## Key Issues Addressed

### 1. emphf Library Segmentation Faults on Ubuntu
**Problem**: The `compute_mphf_seq` binary from the emphf library works correctly on macOS M1 but could cause segmentation faults on Ubuntu x86_64 when built with modified compiler flags.

**Root Cause**: 
- Interference with original build configuration
- Platform-specific compiler optimizations
- SSE instruction differences between ARM64 and x86_64

**Solution Implemented**: 
**Minimal Intervention Approach** - Use the original emphf build process without modifications.

#### Build Process
1. Clone the original emphf repository
2. Navigate to the repository directory  
3. Run `cmake .` (using original CMakeLists.txt)
4. Run `make` (using original build settings)
5. Copy the resulting binaries to our bin directory

This approach ensures:
- No interference with upstream build configuration
- Platform-specific optimizations work as intended
- SSE instructions are handled correctly by the original CMakeLists.txt
- Maximum compatibility across different architectures

## Usage

### Standard Build
```bash
make clean
make all
```

### Cross-Platform Testing
```bash
make test-cross-platform
```

### Safe Build (for problematic systems)
```bash
make clean
make external-safe  # Builds with POPCOUNT disabled
make pybind11       # Continue with rest of build
```

### Debug Platform Issues
```bash
make debug-platform
```

## Build Targets

- `external` - Standard external dependency build using original emphf configuration
- `external-safe` - Conservative build with POPCOUNT disabled for compatibility
- `test-emphf-binary` - Test emphf binary functionality
- `test-cross-platform` - Full cross-platform compatibility test
- `debug-platform` - Display platform and build environment info

## Platform-Specific Notes

### macOS M1 (ARM64)
- Uses original emphf CMakeLists.txt configuration
- Clang compiler with appropriate ARM64 settings
- SSE instructions automatically disabled (not available on ARM)

### Ubuntu x86_64
- Uses original emphf CMakeLists.txt configuration
- GCC compiler with appropriate x86_64 settings
- SSE instructions available and configured by upstream

### Cross-Platform Considerations
- Original emphf build system handles platform differences
- No custom compiler flags interfere with upstream configuration
- Platform detection is handled by CMake in emphf repository

## Troubleshooting

### If emphf binary fails:
1. Try the safe build: `make external-safe`
2. Check platform info: `make debug-platform`
3. Verify binary works: `make test-emphf-binary`

### If build fails entirely:
1. Check CMake version: `cmake --version`
2. Verify compiler: `g++ --version` or `clang++ --version`
3. Check Python environment: `make debug-vars`

## Technical Details

### Why This Approach Works
- **No Build Interference**: We don't modify emphf's CMakeLists.txt or pass custom flags
- **Upstream Maintenance**: Platform-specific optimizations are maintained by emphf developers
- **Architecture Handling**: CMake automatically detects and configures for the target architecture
- **Compiler Compatibility**: Each platform uses its optimal compiler and settings

### Comparison with Previous Approach
- **Before**: Custom compiler flags, architecture detection, manual SSE handling
- **After**: Use original build system, let upstream handle platform differences
- **Result**: Better compatibility, less maintenance, fewer platform-specific bugs

This minimal intervention approach ensures maximum compatibility while requiring minimal maintenance.
