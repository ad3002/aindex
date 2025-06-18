//
// Created by Aleksey Komissarov on 22/01/16.
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <thread>
#include "emphf/common.hpp"
#include <cstdint>
#include <cstring>
#include <limits.h>
#include <math.h>
#include <mutex>
#include <algorithm>
#include <assert.h>
#include "kmers.hpp"
#include <atomic>
#include "hash.hpp"
#include <cassert>
#include "read.hpp"
#include "emphf/common.hpp"
#include <cassert>

int main(int argc, char** argv) {

    if (argc < 8) {
        std::cerr << "Compute AIndex index for genome with pf." << std::endl;
        std::cerr << "Expected arguments: " << argv[0]
        << " <reads_file> <hash_file> <output_prefix> <num_threads> <k> <tf_file> <kmers_bin_file> <kmers_text_file> [index_bin] [indices_bin]" << std::endl;
        std::cerr << "Where:" << std::endl;
        std::cerr << "  reads_file:      Input reads file" << std::endl;
        std::cerr << "  hash_file:       Precomputed emphf hash file" << std::endl;
        std::cerr << "  output_prefix:   Prefix for output files" << std::endl;
        std::cerr << "  num_threads:     Number of threads to use" << std::endl;
        std::cerr << "  k:               K-mer size" << std::endl;
        std::cerr << "  tf_file:         TF file" << std::endl;
        std::cerr << "  kmers_bin_file:  Binary k-mers file" << std::endl;
        std::cerr << "  kmers_text_file: Text k-mers file" << std::endl;
        std::cerr << "  index_bin:       Optional: Output index.bin filename (default: <output_prefix>.index.bin)" << std::endl;
        std::cerr << "  indices_bin:     Optional: Output indices.bin filename (default: <output_prefix>.indices.bin)" << std::endl;
        std::terminate();
    }

    std::vector<READS::READ *> reads;
    PHASH_MAP hash_map;

    std::string read_file = argv[1];
    std::string hash_filename = argv[2]; // precomputed emphf hash
    std::string output_prefix = argv[3];
    static const uint num_threads = atoi(argv[4]);
    Settings::K = atoi(argv[5]);

    std::string tf_file = argv[6];
    std::string kmers_bin_file = argv[7];
    std::string kmers_text_file = argv[8];

    // Optional output filenames
    std::string index_bin_file = (argc > 10) ? argv[10] : output_prefix + ".index.bin";
    std::string indices_bin_file = (argc > 11) ? argv[11] : output_prefix + ".indices.bin";

    emphf::logger() << "Loading hash..." << std::endl;

    if (Settings::K == 13) {
        load_hash_full_tf(hash_map, tf_file, hash_filename);
    } else {
        load_hash(hash_map, hash_filename, tf_file, kmers_bin_file, kmers_text_file);
    }

    emphf::logger() << "\tDone. Kmers: " << hash_map.n << std::endl;

    emphf::logger() << "Load and reads and build docid index..." << std::endl;
    emphf::logger() << "Opening read_file: " << read_file << std::endl;

    std::vector<uint64_t> start_positions;
    std::unordered_map<uint64_t, uint32_t> start2rid;

    std::ifstream infile(read_file);
    if (!infile) {
        emphf::logger() << "Failed open read_file: " << read_file << std::endl;
        exit(10);
    }

    infile.seekg(0, std::ios::end);
    uint64_t length = infile.tellg();
    char *contents = new char[length + 1];
    if (contents == nullptr) {
        emphf::logger() << "Failed to allocate for reads: " << length + 1 << std::endl;
        exit(10);
    }

    std::cout << "Reading reads file: " << read_file;
    infile.seekg(0, std::ios::beg);
    infile.read(contents, length);
    infile.close();
    contents[length] = 0;
    std::cout << " with length: " << length << " bp" << std::endl;
    std::cout << "Done." << std::endl;

    std::cout << "Init aindex..." << std::endl;
    AIndexCompressed aindex(hash_map);
    std::cout << "Done." << std::endl;

    aindex.fill_index_from_reads(contents, length, num_threads, hash_map);

    aindex.save(index_bin_file, indices_bin_file, hash_map);

    emphf::logger() << "\tDone." << std::endl;
    delete[] contents;


    return 0;
}
