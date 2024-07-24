CXX = g++
CXXFLAGS = -std=c++17 -pthread -O3 -fPIC -Wall -Wextra
SRC_DIR = src
INCLUDES = $(SRC_DIR)/helpers.hpp $(SRC_DIR)/debrujin.hpp $(SRC_DIR)/read.hpp $(SRC_DIR)/kmers.hpp $(SRC_DIR)/settings.hpp $(SRC_DIR)/hash.hpp
SOURCES = $(SRC_DIR)/helpers.cpp $(SRC_DIR)/debrujin.cpp $(SRC_DIR)/read.cpp $(SRC_DIR)/kmers.cpp $(SRC_DIR)/settings.cpp $(SRC_DIR)/hash.cpp
OBJECTS = $(SOURCES:.cpp=.o)
BIN_DIR = bin
PACKAGE_DIR = aindex/core
PREFIX = $(CONDA_PREFIX)
INSTALL_DIR = $(PREFIX)/bin

all: clean external $(BIN_DIR) $(BIN_DIR)/compute_index.exe $(BIN_DIR)/compute_aindex.exe $(BIN_DIR)/compute_reads.exe

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/compute_index.exe: $(SRC_DIR)/Compute_index.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/compute_aindex.exe: $(SRC_DIR)/Compute_aindex.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/compute_reads.exe: $(SRC_DIR)/Compute_reads.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(SRC_DIR)/python_wrapper.o: $(SRC_DIR)/python_wrapper.cpp $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(PACKAGE_DIR)/python_wrapper.so: $(SRC_DIR)/python_wrapper.o $(OBJECTS) | $(PACKAGE_DIR)
	$(CXX) $(CXXFLAGS) -shared -Wl,-soname,python_wrapper.so -o $@ $^

external:
	mkdir -p external
	cd external && git clone https://github.com/ad3002/emphf.git
	cd external/emphf && cmake .
	cd external/emphf && make
	cp external/emphf/compute_mphf_seq $(BIN_DIR)
	cp scripts/*py $(BIN_DIR)

install: all
	mkdir -p $(PACKAGE_DIR)
	mkdir -p $(INSTALL_DIR)
	cp bin/compute_index.exe $(INSTALL_DIR)
	cp bin/compute_aindex.exe $(INSTALL_DIR)
	cp bin/compute_reads.exe $(INSTALL_DIR)

clean:
	rm -f $(OBJECTS) $(SRC_DIR)/*.so $(SRC_DIR)/*.o $(BIN_DIR)/*.exe $(PACKAGE_DIR)/python_wrapper.so
	rm -rf external

.PHONY: all clean external install
