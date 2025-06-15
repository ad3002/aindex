# cibuildwheel Configuration Optimization

## Overview

This document describes the optimization of our cibuildwheel configuration based on the official best practices documentation.

## Key Improvements Made

### 1. Configuration Consolidation
- **Before**: Configuration split between `build_wheels.yml` and `pyproject.toml`
- **After**: All cibuildwheel configuration moved to `pyproject.toml` for consistency
- **Benefit**: Single source of truth, easier maintenance

### 2. Simplified CI/CD Workflow
- **Before**: Complex environment variables and platform-specific logic in YAML
- **After**: Clean, minimal `build_wheels.yml` that delegates to `pyproject.toml`
- **Benefit**: Reduced complexity, better maintainability

### 3. Proper Cross-compilation Support
- **Before**: Used `uname -m` which doesn't work for cross-compilation
- **After**: Uses wheel filename analysis to determine target architecture
- **Benefit**: Correct ARM64/x86_64 builds on all runner types

### 4. Platform-specific Optimizations

#### Linux
- Uses package manager detection (`yum` or `apt-get`)
- Builds with `auditwheel` for maximum compatibility
- Supports both x86_64 and aarch64

#### macOS
- Proper MACOSX_DEPLOYMENT_TARGET=10.15 for std::filesystem support
- Architecture detection from wheel filename
- Simplified repair command (copy instead of delocate) for cross-arch compatibility
- Supports both x86_64 and ARM64 with architecture-specific optimizations

#### Windows
- Uses setup.py for builds (more reliable than make on Windows)
- Proper directory creation with Windows syntax
- No external dependencies required

### 5. Build Dependencies Cleanup
- **Before**: Required cmake, git, external emphf repository
- **After**: Only requires standard build tools (make, compiler, Python)
- **Benefit**: Faster builds, no network dependencies, more reliable

## Configuration Structure

### pyproject.toml

```toml
[tool.cibuildwheel]
# Global settings
build = "cp38-* cp39-* cp310-* cp311-* cp312-*"
skip = "*-musllinux*"
test-command = "python -c 'import aindex; print(\"✓ aindex import successful\")'"

# Platform-specific sections
[tool.cibuildwheel.linux]
[tool.cibuildwheel.macos] 
[tool.cibuildwheel.windows]
```

### build_wheels.yml

Simplified to just:
1. Set up Python
2. Install cibuildwheel
3. Run cibuildwheel with config from pyproject.toml
4. Upload artifacts

## Architecture Detection Logic

### macOS Cross-compilation
```bash
WHEEL_FILE='{wheel}' && if echo "$WHEEL_FILE" | grep -q "arm64"; then 
    echo "Building ARM64 version for wheel: $WHEEL_FILE" && make arm64
else 
    echo "Building x86_64 version for wheel: $WHEEL_FILE" && make all
fi
```

This works because cibuildwheel passes the final wheel path to commands, allowing us to detect the target architecture reliably.

## Benefits of the New Configuration

1. **Faster builds**: No external dependencies to download/compile
2. **More reliable**: No network dependencies during build
3. **Better cross-platform support**: Works correctly on all GitHub Actions runners
4. **Easier maintenance**: Single configuration file
5. **Following best practices**: Aligned with cibuildwheel recommendations
6. **Better error handling**: Cleaner error messages and debugging

## Testing

### Local Testing
```bash
# Install cibuildwheel
pip install cibuildwheel

# Test configuration syntax
python -c "import toml; toml.load('pyproject.toml')"

# Test builds work
make clean && make arm64  # or make all
python -c "import aindex; print('✓ Import successful')"
```

### CI Testing
The configuration automatically runs on:
- Pull requests (all platforms)
- Releases (build + publish to PyPI)
- Manual workflow dispatch

## Compatibility

- **Python**: 3.8+ on all platforms
- **Platforms**: Linux (x86_64, aarch64), macOS (x86_64, arm64), Windows (AMD64)
- **Architectures**: Native and cross-compilation supported
- **Dependencies**: Only standard build tools required

## Migration Notes

If updating from previous configuration:
1. Remove cmake/git dependencies
2. Update any custom build scripts to work with new structure
3. Test on target platforms
4. Verify wheel contents and imports work correctly

## References

- [cibuildwheel documentation](https://cibuildwheel.readthedocs.io/)
- [Platform-specific guidance](https://cibuildwheel.readthedocs.io/en/stable/options/)
- [Cross-compilation on macOS](https://cibuildwheel.readthedocs.io/en/stable/options/#macos)
