# Local Build Instructions

This document describes how to build aindex using integrated emphf sources.

## Standard Build

The project includes all necessary emphf sources in `src/emphf/` directory, allowing you to build without external Git dependencies:

```bash
make all
```

This command will:
1. Clean previous builds
2. Copy Python scripts to bin directory
3. Compile all C++ binaries including `compute_mphf_seq` from local sources
4. Build the Python extension
5. Copy all binaries to the package directory

## Key Advantages

- **No external dependencies**: No need to clone the emphf repository
- **Faster builds**: No Git clone or external CMake build steps
- **Cross-platform compatibility**: All sources are integrated and tested
- **ARM64 optimization**: Full support for Apple Silicon with local sources

## Build Targets

- `make all` - Complete build using local emphf sources (DEFAULT)
- `make arm64` - ARM64-optimized build for Apple Silicon
- `make help` - Show all available targets

## Legacy Targets (Deprecated)

- `make all-external` - Build with external emphf repository
- `make simple-all` - Safe build for problematic platforms

## Integrated Files

The following emphf files are integrated in `src/emphf/`:

- `compute_mphf_seq.cpp` - Main target file
- `compute_mphf_generic.hpp` - Core MPHF algorithm
- `hypergraph.hpp` - Hypergraph data structures
- `hypergraph_sorter_seq.hpp` - Sequential sorter
- `mphf.hpp` - Minimal perfect hash function
- `common.hpp` - Common utilities
- `base_hash.hpp` - Hash functions
- And other supporting headers

## Migration Notes

The default `make all` now uses local sources. External dependencies are no longer required.
