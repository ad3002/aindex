# Platform Support

## Supported Platforms (Pre-built Wheels)

AIndex provides pre-built wheels for the following platforms:

### ✅ Primary Support
- **macOS arm64** (Apple Silicon M1/M2/M3)
  - Python 3.8, 3.9, 3.10, 3.11, 3.12
  - Native ARM64 optimizations
  - Full C++ backend functionality
  - Installation: `pip install aindex2`

- **Linux x86_64** 
  - Python 3.8, 3.9, 3.10, 3.11, 3.12
  - Standard Intel/AMD processors
  - Full C++ backend functionality
  - Installation: `pip install aindex2`

## Alternative Installation Methods

### 🔨 Build from Source
For platforms not covered by pre-built wheels:

```bash
git clone https://github.com/ad3002/aindex.git
cd aindex
make all  # or 'make arm64' for Apple Silicon
pip install .
```

**Requirements for building from source:**
- C++ compiler (g++ or clang)
- make
- Python development headers
- pybind11 (automatically installed)

### 🐳 Docker
Use Docker for consistent environment across platforms:

```bash
docker run -it --rm python:3.11-slim bash
pip install aindex2
```

### 💻 Windows Support
Windows users have several options:

1. **Windows Subsystem for Linux (WSL)** - Recommended
   ```bash
   # In WSL Ubuntu terminal
   pip install aindex2
   ```

2. **Docker Desktop**
   ```bash
   docker run -it --rm python:3.11 bash
   pip install aindex2
   ```

3. **Build from source** (requires Visual Studio Build Tools)

### ☁️ Cloud Environments
AIndex works in cloud environments:

- **Google Colab**: `!pip install aindex2`
- **Jupyter Hub/Lab**: `!pip install aindex2`
- **GitHub Codespaces**: `pip install aindex2`

## Performance Notes

- **Apple Silicon (M1/M2/M3)**: Best performance with native ARM64 optimizations
- **Linux x86_64**: Excellent performance with standard optimizations
- **Other platforms**: May have reduced performance when building from source

## Getting Help

If you encounter issues with platform support:

1. Check our [GitHub Issues](https://github.com/ad3002/aindex/issues)
2. Try building from source
3. Use Docker as a fallback option
4. Contact us for specific platform requests

## Future Platform Support

We're continuously working to expand platform support based on user demand. If you need support for a specific platform, please open an issue on GitHub.
