//
// Builds perfect hash for 13-mers - only creates hash file
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <thread>
#include <mutex>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <climits>
#include <cstdlib>
#include <ctime>
#include <sys/stat.h>
#include <errno.h>

#include "emphf/common.hpp"

static std::mutex barrier;

/**
 * Creates perfect hash for 13-mers from k-mers file
 */
void create_13mer_perfect_hash(const std::string& kmers_file, 
                              const std::string& output_hash_file) {
    
    barrier.lock();
    emphf::logger() << "Creating perfect hash for 13-mers..." << std::endl;
    barrier.unlock();
    
    // Verify input file exists
    std::ifstream data_file(kmers_file);
    if (!data_file.is_open()) {
        throw std::runtime_error("Cannot open k-mers file: " + kmers_file);
    }
    
    // Count valid k-mers
    std::string line;
    uint64_t n = 0;
    
    while (std::getline(data_file, line)) {
        if (!line.empty() && line.length() >= 13) {
            std::string kmer = line.substr(0, 13);
            // Check that k-mer contains only ATGC
            bool valid = true;
            for (char c : kmer) {
                if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
                    valid = false;
                    break;
                }
            }
            if (valid) {
                n++;
            }
        }
    }
    data_file.close();
    
    barrier.lock();
    emphf::logger() << "Found " << n << " valid k-mers" << std::endl;
    barrier.unlock();
    
    // Check and create output directory if needed (using POSIX-compatible approach)
    std::string output_path_str = output_hash_file;
    size_t last_slash = output_path_str.find_last_of('/');
    if (last_slash != std::string::npos) {
        std::string output_dir = output_path_str.substr(0, last_slash);
        
        // Check if directory exists using stat
        struct stat st;
        if (stat(output_dir.c_str(), &st) != 0) {
            emphf::logger() << "Output directory does not exist. Creating: " << output_dir << std::endl;
            
            // Create directory recursively
            std::string dir_path = "";
            std::stringstream ss(output_dir);
            std::string segment;
            
            while (std::getline(ss, segment, '/')) {
                if (!segment.empty()) {
                    dir_path += "/" + segment;
                    if (stat(dir_path.c_str(), &st) != 0) {
                        if (mkdir(dir_path.c_str(), 0755) != 0 && errno != EEXIST) {
                            std::cerr << "Error creating directory " << dir_path << ": " << strerror(errno) << std::endl;
                            std::terminate();
                        }
                    }
                }
            }
            emphf::logger() << "Successfully created directory: " << output_dir << std::endl;
        }
    }
    
    // Check if compute_mphf_seq binary exists
    std::ifstream binary_test("./bin/compute_mphf_seq");
    if (!binary_test.is_open()) {
        throw std::runtime_error("Cannot find ./bin/compute_mphf_seq binary");
    }
    binary_test.close();
    
    barrier.lock();
    emphf::logger() << "Building EMPHF hash function using compute_mphf_seq..." << std::endl;
    emphf::logger() << "Processing " << kmers_file << std::endl;
    emphf::logger() << "Memory estimate: ~" << (n * 50 / 1024 / 1024) << " MB needed" << std::endl;
    barrier.unlock();
    
    std::string emphf_cmd = "./bin/compute_mphf_seq " + kmers_file + " " + output_hash_file + " 2>&1";
    barrier.lock();
    emphf::logger() << "Executing: " << emphf_cmd << std::endl;
    barrier.unlock();
    
    int result = system(emphf_cmd.c_str());
    if (result != 0) {
        barrier.lock();
        emphf::logger() << "EMPHF command failed with exit code: " << result << std::endl;
        emphf::logger() << "Checking if output file was created..." << std::endl;
        barrier.unlock();
        
        std::ifstream hash_check(output_hash_file);
        if (!hash_check.is_open()) {
            throw std::runtime_error("EMPHF construction failed - no output file created. This might be a memory issue with " + std::to_string(n) + " k-mers.");
        }
        hash_check.close();
    }
    
    barrier.lock();
    emphf::logger() << "EMPHF hash built successfully, saved to " << output_hash_file << std::endl;
    barrier.unlock();
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <kmers_file> <output_hash_file>" << std::endl;
        std::cerr << "  kmers_file      - File with sequences (one per line)" << std::endl;
        std::cerr << "  output_hash_file - Output hash file (full path)" << std::endl;
        std::cerr << std::endl;
        std::cerr << "Example:" << std::endl;
        std::cerr << "  " << argv[0] << " all_13mers.txt 13mer_index.hash" << std::endl;
        return 1;
    }
    
    std::string kmers_file = argv[1];
    std::string output_hash_file = argv[2];
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "=== 13-mer Perfect Hash Builder ===" << std::endl;
    std::cout << "K-mers file: " << kmers_file << std::endl;
    std::cout << "Output hash file: " << output_hash_file << std::endl;
    std::cout << std::endl;
    
    try {
        create_13mer_perfect_hash(kmers_file, output_hash_file);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Perfect hash creation completed in " << duration.count() << " ms" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}