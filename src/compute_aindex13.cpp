//
// Compute AIndex for k=13 with optimized 13-mer support
// Specialized version for 13-mers using precomputed perfect hash
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

// Specialized structure for 13-mer AIndex
struct AIndex13 {
    uint64_t* indices;        // position indices
    std::atomic<uint64_t>* ppositions; // position completeness
    std::atomic<uint64_t>* positions;  // position itself
    uint64_t total_size = 0;
    uint64_t max_tf = 0;
    static const uint32_t TOTAL_13MERS = 67108864; // 4^13

    AIndex13(const std::string& tf_file) {
        emphf::logger() << "Loading 13-mer frequencies and building index structure..." << std::endl;
        
        // Load tf values
        std::ifstream tf_in(tf_file, std::ios::binary);
        if (!tf_in) {
            std::cerr << "Failed to open tf file: " << tf_file << std::endl;
            exit(1);
        }
        
        std::vector<uint32_t> tf_values(TOTAL_13MERS);
        tf_in.read(reinterpret_cast<char*>(tf_values.data()), TOTAL_13MERS * sizeof(uint32_t));
        tf_in.close();
        
        // Allocate indices array
        emphf::logger() << "Allocating indices array..." << std::endl;
        indices = new uint64_t[TOTAL_13MERS + 1];
        if (indices == nullptr) {
            std::cerr << "Failed to allocate memory for indices: " << TOTAL_13MERS + 1 << std::endl;
            exit(1);
        }
        
        // Build cumulative indices
        indices[0] = 0;
        for (uint64_t i = 1; i < TOTAL_13MERS + 1; ++i) {
            indices[i] = indices[i-1] + tf_values[i-1];
            total_size += tf_values[i-1];
            max_tf = std::max(max_tf, (uint64_t)tf_values[i-1]);
        }
        
        std::cout << "\tmax_tf: " << max_tf << std::endl;
        std::cout << "\ttotal_size: " << total_size << std::endl;
        
        // Allocate position completion counters
        emphf::logger() << "Allocating position completion counters..." << std::endl;
        ppositions = new std::atomic<uint64_t>[TOTAL_13MERS]();
        if (ppositions == nullptr) {
            std::cerr << "Failed to allocate memory for ppositions: " << TOTAL_13MERS << std::endl;
            exit(1);
        }
        
        // Allocate positions array
        emphf::logger() << "Allocating positions array..." << std::endl;
        positions = new std::atomic<uint64_t>[total_size]();
        if (positions == nullptr) {
            std::cerr << "Failed to allocate memory for positions: " << total_size << std::endl;
            exit(1);
        }
        
        emphf::logger() << "AIndex13 structure initialized." << std::endl;
    }
    
    ~AIndex13() {
        if (indices != nullptr) delete[] indices;
        if (ppositions != nullptr) delete[] ppositions;
        if (positions != nullptr) delete[] positions;
    }
    
    void fill_index_from_reads(char* contents, uint64_t length, uint32_t num_threads, 
                              const HASHER& hasher) {
        emphf::logger() << "Building 13-mer index from reads..." << std::endl;
        emphf::logger() << "Length: " << length << ", Threads: " << num_threads << std::endl;
        
        uint64_t batch_size = (length / num_threads) + 1;
        std::vector<std::thread> threads;
        
        for (uint32_t worker_id = 0; worker_id < num_threads; ++worker_id) {
            uint64_t start = worker_id * batch_size;
            uint64_t end = (worker_id + 1) * batch_size;
            if (end > length) {
                end = length;
            }
            
            // Adjust start to avoid cutting k-mers
            if (start > 13) {
                start -= 12;
            }
            
            threads.emplace_back(&AIndex13::worker_thread, this, worker_id, start, end, 
                               contents, std::ref(hasher));
        }
        
        for (auto& t : threads) {
            t.join();
        }
        
        emphf::logger() << "Index building completed." << std::endl;
    }
    
    void worker_thread(int worker_id, uint64_t start, uint64_t end, char* contents, 
                      const HASHER& hasher) {
        emphf::stl_string_adaptor str_adapter;
        static std::mutex progress_mutex;
        
        const int k = 13;
        char ckmer[k + 1];
        
        uint64_t original_start = start;
        uint64_t original_end = end;
        
        // Skip to valid start position
        while (start < end - k + 1) {
            bool found = false;
            for (uint64_t i = start; i < start + k; ++i) {
                if (contents[i] == '\n' || contents[i] == '~' || contents[i] == '?') {
                    start = i + 1;
                    found = true;
                    break;
                }
            }
            if (!found) {
                break;
            }
        }
        
        // Calculate total after adjusting start
        uint64_t total = (end > start + k - 1) ? (end - start - k + 1) : 0;
        int progress = 0;
        int last_progress = -1;  // Initialize to -1 to show 0% immediately
        uint64_t processed = 0;
        uint64_t valid_kmers = 0;
        
        progress_mutex.lock();
        emphf::logger() << "Worker " << worker_id << " started, range: " << original_start 
                       << "-" << original_end << ", adjusted: " << start << "-" << end 
                       << ", total k-mers: " << total << std::endl;
        progress_mutex.unlock();
        
        for (uint64_t i = start; i < end - k + 1; ++i) {
            processed++;
            
            // Safeguard against division by zero and overflow
            if (total > 0) {
                progress = int((100.0 * processed) / total);
                if (progress > 100) progress = 100;  // Cap at 100%
            } else {
                progress = 100;
            }
            
            if (processed % 10000 == 0 || (progress % 10 == 0 && progress != last_progress)) {
                progress_mutex.lock();
                emphf::logger() << "Worker " << worker_id << " progress: " << progress 
                               << "%, processed: " << processed << "/" << total 
                               << ", valid k-mers: " << valid_kmers << std::endl;
                progress_mutex.unlock();
                last_progress = progress;
            }
            
            // Check for invalid characters
            bool skip = false;
            for (int j = 0; j < k; ++j) {
                char c = contents[i + j];
                if (c == '\n' || c == '~' || c == '?' || c == 'N' || 
                    (c != 'A' && c != 'T' && c != 'G' && c != 'C')) {
                    skip = true;
                    break;
                }
            }
            
            if (skip) {
                continue;
            }
            
            // Extract k-mer
            std::memcpy(ckmer, &contents[i], k);
            ckmer[k] = '\0';
            std::string kmer(ckmer);
            
            // Get hash and store position
            try {
                uint64_t hash_id = hasher.lookup(kmer, str_adapter);
                if (hash_id < TOTAL_13MERS) {
                    uint64_t pos_idx = ppositions[hash_id].fetch_add(1, std::memory_order_seq_cst);
                    uint64_t array_idx = indices[hash_id] + pos_idx;
                    
                    if (array_idx < total_size && pos_idx < (indices[hash_id + 1] - indices[hash_id])) {
                        positions[array_idx] = i + 1; // 1-based position
                        valid_kmers++;
                    }
                }
            } catch (...) {
                // Hash lookup failed - try reverse complement
                std::string rev_kmer = get_reverse_complement_13mer(kmer);
                try {
                    uint64_t hash_id = hasher.lookup(rev_kmer, str_adapter);
                    if (hash_id < TOTAL_13MERS) {
                        uint64_t pos_idx = ppositions[hash_id].fetch_add(1, std::memory_order_seq_cst);
                        uint64_t array_idx = indices[hash_id] + pos_idx;
                        
                        if (array_idx < total_size && pos_idx < (indices[hash_id + 1] - indices[hash_id])) {
                            positions[array_idx] = i + 1;
                            valid_kmers++;
                        }
                    }
                } catch (...) {
                    // Both failed - skip
                }
            }
        }
        
        progress_mutex.lock();
        emphf::logger() << "Worker " << worker_id << " finished. Valid k-mers: " << valid_kmers << std::endl;
        progress_mutex.unlock();
    }
    
    std::string get_reverse_complement_13mer(const std::string& kmer) {
        std::string rc = kmer;
        std::reverse(rc.begin(), rc.end());
        for (char& c : rc) {
            switch (c) {
                case 'A': c = 'T'; break;
                case 'T': c = 'A'; break;
                case 'G': c = 'C'; break;
                case 'C': c = 'G'; break;
            }
        }
        return rc;
    }
    
    void get_positions(const std::string& kmer, const HASHER& hasher, std::vector<uint64_t>& result) {
        result.clear();
        
        emphf::stl_string_adaptor str_adapter;
        
        try {
            // Try forward direction
            uint64_t hash_id = hasher.lookup(kmer, str_adapter);
            if (hash_id < TOTAL_13MERS) {
                uint64_t start_idx = indices[hash_id];
                uint64_t end_idx = indices[hash_id + 1];
                
                for (uint64_t i = start_idx; i < end_idx; ++i) {
                    uint64_t pos = positions[i].load();
                    if (pos > 0) { // Valid position
                        result.push_back(pos);
                    }
                }
                return;
            }
        } catch (...) {
            // Try reverse complement
            std::string rev_kmer = get_reverse_complement_13mer(kmer);
            try {
                uint64_t hash_id = hasher.lookup(rev_kmer, str_adapter);
                if (hash_id < TOTAL_13MERS) {
                    uint64_t start_idx = indices[hash_id];
                    uint64_t end_idx = indices[hash_id + 1];
                    
                    for (uint64_t i = start_idx; i < end_idx; ++i) {
                        uint64_t pos = positions[i].load();
                        if (pos > 0) {
                            result.push_back(pos);
                        }
                    }
                }
            } catch (...) {
                // Both failed - return empty
            }
        }
    }
    
    void save(const std::string& index_bin_file, 
              const std::string& indices_bin_file) {
        
        
        
        emphf::logger() << "Saving index.bin array..." << std::endl;
        std::ofstream index_out(index_bin_file, std::ios::binary);
        emphf::logger() << "Positions array size: " << sizeof(uint64_t) * total_size << std::endl;
        
        // Convert atomic to regular array for saving
        std::vector<uint64_t> pos_array(total_size);
        for (uint64_t i = 0; i < total_size; ++i) {
            pos_array[i] = positions[i].load();
        }
        
        index_out.write(reinterpret_cast<const char*>(pos_array.data()), 
                       sizeof(uint64_t) * total_size);
        index_out.close();
        
        emphf::logger() << "Saving indices.bin array..." << std::endl;
        std::ofstream indices_out(indices_bin_file, std::ios::binary);
        indices_out.write(reinterpret_cast<const char*>(indices), 
                         sizeof(uint64_t) * (TOTAL_13MERS + 1));
        indices_out.close();
        
        emphf::logger() << "All files saved successfully." << std::endl;
    }
};

int main(int argc, char** argv) {
    if (argc < 6) {
        std::cerr << "Compute AIndex for 13-mers with perfect hash." << std::endl;
        std::cerr << "Expected arguments: " << argv[0]
                 << " <reads_file> <hash_file> <tf_file> <output_prefix> <num_threads> [pos_bin] [index_bin] [indices_bin]" << std::endl;
        std::cerr << "Where:" << std::endl;
        std::cerr << "  reads_file:      Input reads file (one sequence per line)" << std::endl;
        std::cerr << "  hash_file:       13-mer perfect hash file (.pf)" << std::endl;
        std::cerr << "  tf_file:         13-mer frequencies file (.tf.bin)" << std::endl;
        std::cerr << "  output_prefix:   Prefix for output files" << std::endl;
        std::cerr << "  num_threads:     Number of threads to use" << std::endl;
        std::cerr << "  index_bin:       Optional: Output index.bin filename" << std::endl;
        std::cerr << "  indices_bin:     Optional: Output indices.bin filename" << std::endl;
        return 1;
    }

    std::string reads_file = argv[1];
    std::string hash_file = argv[2];
    std::string tf_file = argv[3];
    std::string output_prefix = argv[4];
    uint32_t num_threads = std::atoi(argv[5]);

    // Optional output filenames
    std::string index_bin_file = (argc > 7) ? argv[7] : output_prefix + ".index.bin";
    std::string indices_bin_file = (argc > 8) ? argv[8] : output_prefix + ".indices.bin";

    emphf::logger() << "Loading 13-mer hash..." << std::endl;
    
    // Load hash
    HASHER hasher;
    std::ifstream hash_in(hash_file, std::ios::binary);
    if (!hash_in) {
        std::cerr << "Failed to open hash file: " << hash_file << std::endl;
        return 1;
    }
    hasher.load(hash_in);
    hash_in.close();
    
    emphf::logger() << "Hash loaded successfully." << std::endl;

    emphf::logger() << "Loading reads and building index..." << std::endl;
    
    // Load reads file
    std::ifstream reads_in(reads_file);
    if (!reads_in) {
        std::cerr << "Failed to open reads file: " << reads_file << std::endl;
        return 1;
    }

    reads_in.seekg(0, std::ios::end);
    uint64_t length = reads_in.tellg();
    char* contents = new char[length + 1];
    if (contents == nullptr) {
        std::cerr << "Failed to allocate memory for reads: " << length + 1 << std::endl;
        return 1;
    }

    reads_in.seekg(0, std::ios::beg);
    reads_in.read(contents, length);
    reads_in.close();
    contents[length] = 0;
    
    emphf::logger() << "Loaded reads of " << length << " bp" << std::endl;

    // Initialize AIndex13 and build index
    emphf::logger() << "Initializing AIndex13..." << std::endl;
    AIndex13 aindex(tf_file);
    
    emphf::logger() << "Building index from reads..." << std::endl;
    aindex.fill_index_from_reads(contents, length, num_threads, hasher);
    
    emphf::logger() << "Saving index..." << std::endl;
    aindex.save(index_bin_file, indices_bin_file);

    delete[] contents;
    
    emphf::logger() << "AIndex13 computation completed successfully!" << std::endl;
    return 0;
}
