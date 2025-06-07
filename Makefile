CXX = g++
CXXFLAGS = -std=c++17 -pthread -O3 -fPIC -Wall -Wextra
LDFLAGS = -shared -Wl,--export-dynamic
SRC_DIR = src
INCLUDES = $(SRC_DIR)/helpers.hpp $(SRC_DIR)/debrujin.hpp $(SRC_DIR)/read.hpp $(SRC_DIR)/kmers.hpp $(SRC_DIR)/settings.hpp $(SRC_DIR)/hash.hpp
SOURCES = $(SRC_DIR)/helpers.cpp $(SRC_DIR)/debrujin.cpp $(SRC_DIR)/read.cpp $(SRC_DIR)/kmers.cpp $(SRC_DIR)/settings.cpp $(SRC_DIR)/hash.cpp
OBJECTS = $(SOURCES:.cpp=.o)
BIN_DIR = bin
PACKAGE_DIR = aindex/core
PREFIX = $(CONDA_PREFIX)
INSTALL_DIR = $(PREFIX)/bin

all: clean external $(BIN_DIR) $(BIN_DIR)/compute_index.exe $(BIN_DIR)/compute_aindex.exe $(BIN_DIR)/compute_reads.exe $(PACKAGE_DIR)/python_wrapper.so

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/compute_index.exe: $(SRC_DIR)/Compute_index.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/compute_aindex.exe: $(SRC_DIR)/Compute_aindex.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/compute_reads.exe: $(SRC_DIR)/Compute_reads.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $< -o $@

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
	fi
	@echo "Building emphf..."
	@echo "Checking for cmake..."
	@if command -v /usr/bin/cmake >/dev/null 2>&1; then \
		echo "Using system cmake: /usr/bin/cmake"; \
		cd external/emphf && /usr/bin/cmake . && make; \
	elif python3 -c "import cmake" >/dev/null 2>&1; then \
		echo "Python cmake module found, installing system cmake..."; \
		apt-get update && apt-get install -y cmake; \
		cd external/emphf && /usr/bin/cmake . && make; \
	else \
		echo "Installing cmake..."; \
		apt-get update && apt-get install -y cmake; \
		cd external/emphf && cmake . && make; \
	fi || { \
		echo "Failed to build emphf. Trying alternative approach..."; \
		echo "Removing Python cmake and installing system cmake..."; \
		pip uninstall -y cmake || true; \
		apt-get update && apt-get install -y cmake; \
		cd external/emphf && cmake . && make; \
	}
	cp external/emphf/compute_mphf_seq $(BIN_DIR)/ || { \
		echo "Failed to copy emphf binary. Build may have failed."; \
		exit 1; \
	}
	cp scripts/compute_aindex.py $(BIN_DIR)/
	cp scripts/compute_index.py $(BIN_DIR)/
	cp scripts/reads_to_fasta.py $(BIN_DIR)/
	@echo "External dependencies setup complete."

install: all
	mkdir -p ${BIN_DIR}
	mkdir -p $(PACKAGE_DIR)
	mkdir -p $(INSTALL_DIR)
	cp bin/compute_index.exe $(INSTALL_DIR)/
	cp bin/compute_aindex.exe $(INSTALL_DIR)/
	cp bin/compute_reads.exe $(INSTALL_DIR)/

clean:
	rm -f $(OBJECTS) $(SRC_DIR)/*.so $(SRC_DIR)/*.o $(BIN_DIR)/*.exe $(PACKAGE_DIR)/python_wrapper.so
	rm -rf external

.PHONY: all clean external install