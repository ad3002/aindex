# Local Build Instructions

This document describes how to build aindex using integrated emphf sources instead of external dependencies.

## Local Build (Recommended)

The project now includes all necessary emphf sources in `src/emphf/` directory, allowing you to build without external Git dependencies:

```bash
make all-local
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

- `make all-local` - Complete build using local emphf sources (RECOMMENDED)
- `make all` - Traditional build with external emphf repository
- `make help` - Show all available targets

## Integrated Files

The following emphf files are now integrated in `src/emphf/`:

- `compute_mphf_seq.cpp` - Main target file
- `compute_mphf_generic.hpp` - Core MPHF algorithm
- `hypergraph.hpp` - Hypergraph data structures
- `hypergraph_sorter_seq.hpp` - Sequential sorter
- `mphf.hpp` - Minimal perfect hash function
- `common.hpp` - Common utilities
- `base_hash.hpp` - Hash functions
- And other supporting headers

## Migration Notes

If you were previously using `make all`, you can now use `make all-local` for the same functionality without external dependencies.

The local build produces identical results to the external build but is more reliable and faster.
