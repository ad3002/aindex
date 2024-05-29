# Makefile for aindex project

CXX = g++
CXXFLAGS = -std=c++17 -pthread -O3
SRC_DIR = src
INCLUDES = $(SRC_DIR)/helpers.hpp $(SRC_DIR)/debrujin.hpp $(SRC_DIR)/read.hpp $(SRC_DIR)/kmers.hpp $(SRC_DIR)/settings.hpp $(SRC_DIR)/hash.hpp
SOURCES = $(SRC_DIR)/helpers.cpp $(SRC_DIR)/debrujin.cpp $(SRC_DIR)/read.cpp $(SRC_DIR)/kmers.cpp $(SRC_DIR)/settings.cpp $(SRC_DIR)/hash.cpp
OBJECTS = $(SOURCES:.cpp=.o)
BIN_DIR = bin

all: $(BIN_DIR) $(BIN_DIR)/compute_index.exe $(BIN_DIR)/compute_aindex.exe $(BIN_DIR)/compute_reads.exe python_wrapper_linux.so python_wrapper_mac.so

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/compute_index.exe: $(SRC_DIR)/Compute_index.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/compute_aindex.exe: $(SRC_DIR)/Compute_aindex.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/compute_reads.exe: $(SRC_DIR)/Compute_reads.cpp $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

python_wrapper_linux.so: $(SRC_DIR)/python_wrapper.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) -shared -Wl,-soname,python_wrapper.so -o $@ $^

python_wrapper_mac.so: $(SRC_DIR)/python_wrapper.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) -shared -Wl,-install_name,python_wrapper.so -o $@ $^

$(SRC_DIR)/python_wrapper.o: $(SRC_DIR)/python_wrapper.cpp $(INCLUDES)
	$(CXX) $(CXXFLAGS) -fPIC -c $< -o $@

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(SRC_DIR)/*.o $(BIN_DIR)/*.exe python_wrapper_linux.so python_wrapper_mac.so

.PHONY: all clean
