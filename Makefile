CXX = g++
CXXFLAGS = -std=c++17 -pthread -O3 -fPIC -Wall -Wextra
LDFLAGS = -shared -Wl,--export-dynamic
SRC_DIR = src
OBJ_DIR = obj
INCLUDES = $(SRC_DIR)/helpers.hpp $(SRC_DIR)/debrujin.hpp $(SRC_DIR)/read.hpp $(SRC_DIR)/kmers.hpp $(SRC_DIR)/settings.hpp $(SRC_DIR)/hash.hpp
SOURCES = $(SRC_DIR)/helpers.cpp $(SRC_DIR)/debrujin.cpp $(SRC_DIR)/read.cpp $(SRC_DIR)/kmers.cpp $(SRC_DIR)/settings.cpp $(SRC_DIR)/hash.cpp
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
BIN_DIR = bin
PACKAGE_DIR = aindex/core
PREFIX = $(CONDA_PREFIX)
INSTALL_DIR = $(PREFIX)/bin

# Python and pybind11 configuration - auto-detect available Python
# In cibuildwheel, use the current Python; otherwise try to find the best one
PYTHON_CMD := $(shell \
    if [ -n "$$CIBUILDWHEEL" ] && which python >/dev/null 2>&1; then \
        echo python; \
    elif /opt/homebrew/opt/python@3.11/bin/python3.11 --version >/dev/null 2>&1; then \
        echo /opt/homebrew/opt/python@3.11/bin/python3.11; \
    elif python3.11 --version >/dev/null 2>&1; then \
        echo python3.11; \
    elif python3 --version >/dev/null 2>&1; then \
        echo python3; \
    elif python --version >/dev/null 2>&1; then \
        echo python; \
    else \
        echo python3; \
    fi)

# Auto-detect Python version and corresponding config
PYTHON_VERSION := $(shell $(PYTHON_CMD) -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")

# Try to find the correct python-config for the active Python
# In cibuildwheel, try the most common patterns first
PYTHON_CONFIG_CANDIDATES := $(shell \
    if [ -n "$$CIBUILDWHEEL" ]; then \
        echo "$(PYTHON_CMD)-config python$(PYTHON_VERSION)-config python3-config python-config"; \
    else \
        echo "$(PYTHON_CMD)-config python$(PYTHON_VERSION)-config python3-config python-config"; \
    fi)
PYTHON_CONFIG = $(shell for cmd in $(PYTHON_CONFIG_CANDIDATES); do \
    if which $$cmd >/dev/null 2>&1; then \
        echo $$cmd; \
        break; \
    fi; \
done)

# Fallback if no config found
ifeq ($(PYTHON_CONFIG),)
    PYTHON_CONFIG = python3-config
endif

# Safely get pybind11 include path
PYTHON_INCLUDE := $(shell $(PYTHON_CMD) -c "try: import pybind11; print(pybind11.get_include());\nexcept: print('')" 2>/dev/null)
PYTHON_HEADERS := $(shell $(PYTHON_CONFIG) --includes)
PYTHON_SUFFIX := $(shell $(PYTHON_CONFIG) --extension-suffix)

# Detect OS for macOS-specific settings
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    CXXFLAGS += -stdlib=libc++ -I$(PYTHON_INCLUDE) $(PYTHON_HEADERS)
    LDFLAGS = -shared -undefined dynamic_lookup
    MACOS = true
else
    CXXFLAGS += -I$(PYTHON_INCLUDE) $(PYTHON_HEADERS)
    LDFLAGS = -shared
    MACOS = false
endif

# Platform-specific binary extensions
ifeq ($(UNAME_S),Windows_NT)
    BIN_EXT = .exe
else
    BIN_EXT = 
endif

# Binary targets with platform-appropriate extensions
BINARIES = $(BIN_DIR)/compute_index$(BIN_EXT) $(BIN_DIR)/compute_aindex$(BIN_EXT) $(BIN_DIR)/compute_reads$(BIN_EXT) $(BIN_DIR)/kmer_counter$(BIN_EXT) $(BIN_DIR)/generate_all_13mers$(BIN_EXT) $(BIN_DIR)/build_13mer_hash$(BIN_EXT) $(BIN_DIR)/count_kmers13$(BIN_EXT) $(BIN_DIR)/compute_aindex13$(BIN_EXT)

all: clean external $(BIN_DIR) $(OBJ_DIR) $(BINARIES) pybind11 copy-to-package

# Alternative simplified all target that matches what's used in setup.py
simple-all: clean external-safe $(BIN_DIR) $(OBJ_DIR) $(BINARIES) pybind11 copy-to-package

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR)/compute_index$(BIN_EXT): $(SRC_DIR)/Compute_index.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/compute_aindex$(BIN_EXT): $(SRC_DIR)/Compute_aindex.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/compute_reads$(BIN_EXT): $(SRC_DIR)/Compute_reads.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/kmer_counter$(BIN_EXT): $(SRC_DIR)/Count_kmers.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN_DIR)/generate_all_13mers$(BIN_EXT): $(SRC_DIR)/generate_all_13mers.cpp $(OBJ_DIR)/kmers.o | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/build_13mer_hash$(BIN_EXT): $(SRC_DIR)/Build_13mer_hash.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I./external $^ -o $@

$(BIN_DIR)/count_kmers13$(BIN_EXT): $(SRC_DIR)/count_kmers13.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I./external $^ -o $@

$(BIN_DIR)/compute_aindex13$(BIN_EXT): $(SRC_DIR)/Compute_aindex13.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INC) -I./external $< $(OBJECTS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INCLUDES) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Pybind11 module
pybind11: $(OBJECTS) $(SRC_DIR)/python_wrapper.cpp | $(PACKAGE_DIR)
	@echo "=== Building Python extension ==="
	@echo "Python command: $(PYTHON_CMD)"
	@echo "Python version: $(PYTHON_VERSION)"
	@echo "Python config: $(PYTHON_CONFIG)"
	@echo "Extension suffix: $(PYTHON_SUFFIX)"
	@echo "CIBUILDWHEEL env: $$CIBUILDWHEEL"
	@$(PYTHON_CMD) -c "import sys; print(f'Active Python: {sys.executable}')"
	@PYBIND11_INCLUDE=$$($(PYTHON_CMD) -c "import pybind11; print(pybind11.get_include())" 2>/dev/null) && \
	if [ -z "$$PYBIND11_INCLUDE" ]; then \
		echo "Error: pybind11 not found. Please install pybind11: pip install pybind11"; \
		exit 1; \
	else \
		echo "pybind11 include path: $$PYBIND11_INCLUDE"; \
		$(CXX) $(CXXFLAGS) -I$$PYBIND11_INCLUDE -I./external $(LDFLAGS) -o $(PACKAGE_DIR)/aindex_cpp$(PYTHON_SUFFIX) $(SRC_DIR)/python_wrapper.cpp $(OBJECTS); \
	fi

$(PACKAGE_DIR)/python_wrapper.so: $(SRC_DIR)/python_wrapper.o $(OBJECTS) | $(PACKAGE_DIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

external:
	@echo "Setting up external dependencies..."
	mkdir -p ${BIN_DIR}
	mkdir -p external
	mkdir -p $(PACKAGE_DIR)
	@if [ ! -d "external/emphf" ]; then \
		echo "Cloning emphf repository..."; \
		cd external && git clone https://github.com/ad3002/emphf.git || { \
			echo "Failed to clone emphf repository. Please check your internet connection."; \
			exit 1; \
		}; \
		echo "Applying CMake version patch..."; \
		cd emphf && patch -p1 < ../../patches/emphf_cmake_version.patch; \
	fi
	@echo "Building emphf using original build process..."
	@echo "Platform: $(shell uname -s) $(shell uname -m)"
	cd external/emphf && env -i PATH="$$PATH" HOME="$$HOME" cmake . && env -i PATH="$$PATH" HOME="$$HOME" make
	@echo "Copying emphf binaries to our bin directory..."
	@if [ -f "external/emphf/compute_mphf_seq" ]; then \
		echo "✓ compute_mphf_seq found"; \
		cp external/emphf/compute_mphf_seq $(BIN_DIR)/; \
		echo "✓ Binary copied to $(BIN_DIR)/"; \
	else \
		echo "✗ compute_mphf_seq not found"; \
		ls -la external/emphf/compute_mphf* || echo "No compute_mphf* files found"; \
		exit 1; \
	fi
	@echo "Copying our Python scripts..."
	cp scripts/compute_aindex.py $(BIN_DIR)/
	cp scripts/compute_index.py $(BIN_DIR)/
	cp scripts/reads_to_fasta.py $(BIN_DIR)/
	@echo "External dependencies setup complete."

# Alternative build for problematic platforms
external-safe: 
	@echo "Setting up external dependencies with safe mode..."
	mkdir -p ${BIN_DIR}
	mkdir -p external
	mkdir -p $(PACKAGE_DIR)
	@if [ ! -d "external/emphf" ]; then \
		echo "Cloning emphf repository..."; \
		cd external && git clone https://github.com/ad3002/emphf.git || { \
			echo "Failed to clone emphf repository. Please check your internet connection."; \
			exit 1; \
		}; \
		echo "Applying CMake version patch..."; \
		cd emphf && patch -p1 < ../../patches/emphf_cmake_version.patch; \
	fi
	@echo "Building emphf with POPCOUNT disabled for compatibility..."
	@echo "Platform: $(shell uname -s) $(shell uname -m)"
	cd external/emphf && env -i PATH="$$PATH" HOME="$$HOME" cmake -DEMPHF_USE_POPCOUNT=OFF . && env -i PATH="$$PATH" HOME="$$HOME" make
	@echo "Copying emphf binaries to our bin directory..."
	@if [ -f "external/emphf/compute_mphf_seq" ]; then \
		echo "✓ compute_mphf_seq found"; \
		cp external/emphf/compute_mphf_seq $(BIN_DIR)/; \
		echo "✓ Binary copied to $(BIN_DIR)/"; \
	else \
		echo "✗ compute_mphf_seq not found"; \
		ls -la external/emphf/compute_mphf* || echo "No compute_mphf* files found"; \
		exit 1; \
	fi
	@echo "Copying our Python scripts..."
	cp scripts/compute_aindex.py $(BIN_DIR)/
	cp scripts/compute_index.py $(BIN_DIR)/
	cp scripts/reads_to_fasta.py $(BIN_DIR)/
	@echo "Safe external dependencies setup complete."

install: all
	mkdir -p ${BIN_DIR}
	mkdir -p $(PACKAGE_DIR)
	mkdir -p $(INSTALL_DIR)
	cp bin/compute_index.exe $(INSTALL_DIR)/
	cp bin/compute_aindex.exe $(INSTALL_DIR)/
	cp bin/compute_reads.exe $(INSTALL_DIR)/
	cp bin/kmer_counter.exe $(INSTALL_DIR)/
	cp bin/generate_all_13mers.exe $(INSTALL_DIR)/
	cp bin/Build_13mer_hash.exe $(INSTALL_DIR)/
	cp bin/compute_aindex13.exe $(INSTALL_DIR)/
	cp bin/count_kmers13.exe $(INSTALL_DIR)/

clean:
	rm -rf $(OBJ_DIR) $(SRC_DIR)/*.so $(BIN_DIR) $(PACKAGE_DIR)/python_wrapper.so $(PACKAGE_DIR)/aindex_cpp*.so
	rm -rf external
	rm -rf aindex/bin

# macOS-specific target for manual compilation
macos: clean $(PACKAGE_DIR) $(OBJ_DIR)
	@echo "Building for macOS..."
	@echo "Note: This target skips external dependencies due to ARM64 compatibility issues"
	mkdir -p $(PACKAGE_DIR)
	mkdir -p $(OBJ_DIR)
	g++ -c -std=c++11 -fPIC $(SRC_DIR)/python_wrapper.cpp -o $(OBJ_DIR)/python_wrapper.o
	g++ -c -std=c++11 -fPIC $(SRC_DIR)/kmers.cpp -o $(OBJ_DIR)/kmers.o
	g++ -c -std=c++11 -fPIC $(SRC_DIR)/debrujin.cpp -o $(OBJ_DIR)/debrujin.o
	g++ -c -std=c++11 -fPIC $(SRC_DIR)/hash.cpp -o $(OBJ_DIR)/hash.o
	g++ -c -std=c++11 -fPIC $(SRC_DIR)/read.cpp -o $(OBJ_DIR)/read.o
	g++ -c -std=c++11 -fPIC $(SRC_DIR)/settings.cpp -o $(OBJ_DIR)/settings.o
	g++ -c -std=c++11 -fPIC $(SRC_DIR)/helpers.cpp -o $(OBJ_DIR)/helpers.o
	g++ -shared -Wl,-install_name,python_wrapper.so -o $(PACKAGE_DIR)/python_wrapper.so \
		$(OBJ_DIR)/python_wrapper.o $(OBJ_DIR)/kmers.o $(OBJ_DIR)/debrujin.o $(OBJ_DIR)/hash.o $(OBJ_DIR)/read.o $(OBJ_DIR)/settings.o $(OBJ_DIR)/helpers.o
	@echo "macOS build complete! python_wrapper.so created in $(PACKAGE_DIR)/"

# macOS simplified target for testing without emphf dependencies
macos-simple: clean $(PACKAGE_DIR) $(OBJ_DIR)
	@echo "Building simplified version for macOS (testing only)..."
	mkdir -p $(PACKAGE_DIR)
	mkdir -p $(OBJ_DIR)
	g++ -c -std=c++11 -fPIC $(SRC_DIR)/python_wrapper_simple.cpp -o $(OBJ_DIR)/python_wrapper_simple.o
	g++ -shared -Wl,-install_name,python_wrapper.so -o $(PACKAGE_DIR)/python_wrapper.so \
		$(OBJ_DIR)/python_wrapper_simple.o
	@echo "macOS simplified build complete! python_wrapper.so created in $(PACKAGE_DIR)/"

# Create package directory
$(PACKAGE_DIR):
	mkdir -p $(PACKAGE_DIR)

# Copy binaries to package directory for inclusion in wheel
copy-to-package: $(BINARIES)
	@echo "Copying binaries to package directory..."
	mkdir -p aindex/bin
	cp -f $(BIN_DIR)/* aindex/bin/ 2>/dev/null || true
	@echo "✓ Binaries copied to aindex/bin/"
	@ls -la aindex/bin/ || echo "No files in aindex/bin/"

# Test targets
test:
	@echo "Running full regression tests..."
	$(PYTHON_CMD) test_aindex_functionality.py
	$(PYTHON_CMD) test_aindex_functionality_k13.py

test-all:
	@echo "Running full regression tests..."
	$(PYTHON_CMD) test_aindex_functionality.py
	$(PYTHON_CMD) test_aindex_functionality_k13.py

# Debug target for cross-platform issues
debug-platform:
	@echo "=== Platform Debug Information ==="
	@echo "OS: $(shell uname -s)"
	@echo "Architecture: $(shell uname -m)"
	@echo "Compiler: $(CXX)"
	@echo "C++ Standard Library:"
	@$(CXX) --version || true
	@echo "Python: $(PYTHON_CMD)"
	@$(PYTHON_CMD) --version
	@echo "Python Config: $(PYTHON_CONFIG)"
	@echo "CMake version:"
	@cmake --version 2>/dev/null || echo "CMake not found"
	@echo "=== Build Flags ==="
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
	@echo "=================================="

# Help target
help:
	@echo "Available targets:"
	@echo "  all              - Build all binaries and Python extension"
	@echo "  simple-all       - Build with safe external dependencies"
	@echo "  clean            - Clean build artifacts"
	@echo "  pybind11         - Build only the Python extension"
	@echo "  test             - Run Python API tests"
	@echo "  test-all         - Run Python API tests"
	@echo "  debug-platform   - Display platform and build environment information"
	@echo "  install          - Install binaries to system (requires CONDA_PREFIX)"
	@echo "  help             - Show this help message"
	@echo ""
	@echo "Recommended usage:"
	@echo "  make all         - Complete build (external dependencies + binaries + Python)"
	@echo "  make simple-all  - Safe build for problematic platforms"
	@echo "  make test-all    - Complete test suite for new users/CI"
	@echo ""
	@echo "Cross-platform debugging:"
	@echo "  make debug-platform - Show platform information"
	@echo "  make external-safe - Safe build for problematic platforms"
	@echo ""
	@echo "Documentation:"
	@echo "  See CROSS_PLATFORM.md for cross-platform compatibility details"
	@echo ""
	@echo "Python version detected: $(PYTHON_VERSION)"
	@echo "Python config: $(PYTHON_CONFIG)"
	@echo "Extension suffix: $(PYTHON_SUFFIX)"

# Debug target to print variables
debug-vars:
	@echo "PYTHON_CMD: $(PYTHON_CMD)"
	@echo "PYTHON_VERSION: $(PYTHON_VERSION)"
	@echo "PYTHON_CONFIG: $(PYTHON_CONFIG)"
	@echo "PYTHON_INCLUDE: $(PYTHON_INCLUDE)"
	@echo "PYTHON_HEADERS: $(PYTHON_HEADERS)"
	@echo "PYTHON_SUFFIX: $(PYTHON_SUFFIX)"

.PHONY: all simple-all clean external external-safe install macos macos-simple test test-all debug-platform help debug-vars copy-to-package