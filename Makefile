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
# Try to find the actual Python executable, handle aliases properly
PYTHON_CMD := $(shell \
    if /opt/homebrew/opt/python@3.11/bin/python3.11 --version >/dev/null 2>&1; then \
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
# First try with the exact Python command that's being used
PYTHON_CONFIG_CANDIDATES = $(PYTHON_CMD)-config python$(PYTHON_VERSION)-config python3-config python-config
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

PYTHON_INCLUDE := $(shell $(PYTHON_CMD) -c "import pybind11; print(pybind11.get_include())")
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
    MACOS = false
endif

all: clean external $(BIN_DIR) $(OBJ_DIR) $(BIN_DIR)/compute_index.exe $(BIN_DIR)/compute_aindex.exe $(BIN_DIR)/compute_reads.exe $(BIN_DIR)/kmer_counter.exe $(BIN_DIR)/generate_all_13mers.exe $(BIN_DIR)/Build_13mer_hash.exe $(BIN_DIR)/count_kmers13.exe $(BIN_DIR)/compute_aindex13.exe pybind11

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR)/compute_index.exe: $(SRC_DIR)/Compute_index.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/compute_aindex.exe: $(SRC_DIR)/Compute_aindex.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/compute_reads.exe: $(SRC_DIR)/Compute_reads.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/kmer_counter.exe: $(SRC_DIR)/Count_kmers.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN_DIR)/generate_all_13mers.exe: $(SRC_DIR)/generate_all_13mers.cpp $(OBJ_DIR)/kmers.o | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/Build_13mer_hash.exe: $(SRC_DIR)/Build_13mer_hash.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I./external $^ -o $@

$(BIN_DIR)/count_kmers13.exe: $(SRC_DIR)/count_kmers13.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I./external $^ -o $@

$(BIN_DIR)/compute_aindex13.exe: $(SRC_DIR)/Compute_aindex13.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INC) -I./external $< $(OBJECTS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INCLUDES) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Pybind11 module
pybind11: $(OBJECTS) $(SRC_DIR)/python_wrapper.cpp | $(PACKAGE_DIR)
	@echo "Building Python extension for Python $(PYTHON_VERSION)"
	@echo "Using Python config: $(PYTHON_CONFIG)"
	@echo "Extension suffix: $(PYTHON_SUFFIX)"
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(PACKAGE_DIR)/aindex_cpp$(PYTHON_SUFFIX) $(SRC_DIR)/python_wrapper.cpp $(OBJECTS)

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

# Test targets
test: test-python-api

test-python-api: pybind11
	@echo "Running Python API tests..."
	@$(PYTHON_CMD) tests/test_python_api_basic.py

# Test targets without dependencies (for use in test-all)
test-python-api-only:
	@echo "Running Python API tests..."
	@$(PYTHON_CMD) tests/test_python_api_basic.py

test-python-with-data: pybind11
	@echo "Running Python API integration tests..."
	@$(PYTHON_CMD) tests/test_python_api_integration.py

test-python-with-data-only:
	@echo "Running Python API integration tests..."
	@$(PYTHON_CMD) tests/test_python_api_integration.py

test-regression: all
	@echo "Running full regression tests..."
	@cd tests && $(PYTHON_CMD) tests/test_regression.py --skip-jellyfish

test-regression-only:
	@echo "Running full regression tests..."
	@cd tests && $(PYTHON_CMD) tests/test_regression.py --skip-jellyfish

test-quick: pybind11
	@echo "Running quick validation test..."
	@$(PYTHON_CMD) tests/quick_test.py

test-demo: pybind11 test-regression-only
	@echo "Running demo script with test data..."
	@$(PYTHON_CMD) tests/demo.py

test-demo-comprehensive: pybind11 test-regression-only
	@echo "Running comprehensive demo script with test data..."
	@$(PYTHON_CMD) tests/demo.py --comprehensive

test-demo-only:
	@echo "Running demo script with existing test data..."
	@$(PYTHON_CMD) tests/demo.py

test-demo-comprehensive-only:
	@echo "Running comprehensive demo script with existing test data..."
	@$(PYTHON_CMD) tests/demo.py --comprehensive

test-performance: pybind11 test-regression-only
	@echo "Running performance benchmark tests..."
	@$(PYTHON_CMD) tests/performance_benchmark.py

test-performance-only:
	@echo "Running performance benchmark tests with existing data..."
	@$(PYTHON_CMD) tests/performance_benchmark.py

test-speed: pybind11 test-regression-only
	@echo "Running speed tests..."
	@$(PYTHON_CMD) tests/speed_test.py

test-speed-only:
	@echo "Running speed tests with existing data..."
	@$(PYTHON_CMD) tests/speed_test.py

test-full: test-regression test-python-with-data
	@echo "🎉 All tests completed successfully!"

test-all: all
	@echo "=========================================="
	@echo "🧪 RUNNING ALL AINDEX TESTS"
	@echo "=========================================="
	@echo "Project already built. Running tests..."
	@echo ""
	@echo "Step 1/2: Running regression tests (generates test data)..."
	@$(MAKE) --no-print-directory test-regression-only
	@echo ""
	@echo "Step 2/2: Running Python API tests..."
	@$(MAKE) --no-print-directory test-python-api-only
	@$(MAKE) --no-print-directory test-python-with-data-only
	@echo ""
	@echo "=========================================="
	@echo "🎉 ALL TESTS COMPLETED SUCCESSFULLY! 🎉"
	@echo "=========================================="

# Cross-platform testing targets
test-emphf-binary: external
	@echo "Testing emphf binary compatibility..."
	@echo "Platform: $(shell uname -s) $(shell uname -m)"
	@if [ -x "$(BIN_DIR)/compute_mphf_seq" ]; then \
		echo "✓ Binary is executable"; \
		echo "Testing with minimal input..."; \
		echo -e "test1\ntest2\ntest3" > /tmp/test_input.txt; \
		if $(BIN_DIR)/compute_mphf_seq /tmp/test_input.txt /tmp/test_output.mphf 2>/dev/null; then \
			echo "✓ Binary runs successfully on $(shell uname -s) $(shell uname -m)"; \
			rm -f /tmp/test_input.txt /tmp/test_output.mphf; \
		else \
			echo "✗ Binary fails on $(shell uname -s) $(shell uname -m)"; \
			rm -f /tmp/test_input.txt /tmp/test_output.mphf; \
			exit 1; \
		fi; \
	else \
		echo "✗ Binary not found or not executable"; \
		exit 1; \
	fi

test-cross-platform: test-emphf-binary
	@echo "Cross-platform compatibility test passed!"

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
	@echo "  clean            - Clean build artifacts"
	@echo "  pybind11         - Build only the Python extension"
	@echo "  test             - Run basic Python API tests (alias for test-python-api)"
	@echo "  test-quick       - Quick validation test (fast functionality check)"
	@echo "  test-demo        - Run demo script (builds test data if needed)"
	@echo "  test-demo-comprehensive - Run comprehensive demo (builds test data if needed)"
	@echo "  test-performance - Run comprehensive performance benchmarks"
	@echo "  test-speed       - Run quick speed tests (key metrics)"
	@echo "  test-python-api  - Test Python module import and basic functionality"
	@echo "  test-python-with-data - Test Python API with real data (requires test data)"
	@echo "  test-regression  - Run full regression tests and generate test data"
	@echo "  test-full        - Run all tests (regression + Python API)"
	@echo "  test-all         - Run complete test suite (build once + all tests)"
	@echo "  test-emphf-binary - Test emphf binary compatibility on current platform"
	@echo "  test-cross-platform - Test cross-platform compatibility"
	@echo "  debug-platform   - Display platform and build environment information"
	@echo "  install          - Install binaries to system (requires CONDA_PREFIX)"
	@echo "  help             - Show this help message"
	@echo ""
	@echo "Recommended usage:"
	@echo "  make test-all    - Complete test suite for new users/CI (optimized - no rebuild)"
	@echo "  make test-cross-platform - Test cross-platform compatibility"
	@echo "  make test-demo   - Interactive demo with real data"
	@echo "  make test-demo-comprehensive - Full functionality showcase"
	@echo "  make test-speed  - Quick performance check (key metrics)"
	@echo "  make test-performance - Comprehensive performance analysis"
	@echo "  make test-quick  - Quick validation for development"
	@echo "  make test        - Basic Python API tests for development"
	@echo ""
	@echo "Cross-platform debugging:"
	@echo "  make debug-platform - Show platform information"
	@echo "  make test-emphf-binary - Test emphf binary specifically"
	@echo "  make external-safe - Safe build for problematic platforms"
	@echo ""
	@echo "Documentation:"
	@echo "  See CROSS_PLATFORM.md for cross-platform compatibility details"
	@echo ""
	@echo "Manual demo usage:"
	@echo "  python tests/demo.py                - Basic demo"
	@echo "  python tests/demo.py --comprehensive - Full comprehensive demo"
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

.PHONY: all clean external external-safe install macos macos-simple test test-quick test-demo test-demo-comprehensive test-demo-only test-demo-comprehensive-only test-python-api test-python-api-only test-python-with-data test-python-with-data-only test-regression test-regression-only test-full test-all test-emphf-binary test-cross-platform debug-platform help debug-vars