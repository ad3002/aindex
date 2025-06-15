//
// ARM64-optimized Compute AIndex for k=13 with optimized 13-mer support
// Optimized for Apple Silicon (M1/M2/M3) with reduced memory usage during compilation
// Specialized version for 13-mers using precomputed perfect hash
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <unordered_map>

#ifdef __aarch64__
#include <arm_neon.h>
#endif

// Minimal includes for ARM64 to reduce compilation memory usage
#include "hash.hpp"
#include "read.hpp"

// ARM64-optimized structure for 13-mer AIndex
class ARM64AIndex13 {
private:
    static const uint32_t TOTAL_13MERS = 67108864; // 4^13
    static const size_t CACHE_LINE_SIZE = 128; // ARM64 cache line
    static const size_t BATCH_SIZE = 10000; // Process k-mers in batches
    
    uint64_t* indices;
    std::atomic<uint32_t>* ppositions; // Use 32-bit counters to save memory
    std::atomic<uint64_t>* positions;
    uint64_t total_size = 0;
    uint64_t max_tf = 0;
    
    // ARM64 optimization: lookup tables
    static const bool valid_nucleotide[256];
    static const unsigned char complement_table[256];
    
public:
    ARM64AIndex13(const std::string& tf_file) {
        std::cout << "[ARM64] Loading 13-mer frequencies and building index structure..." << std::endl;
        
        // Load tf values efficiently
        std::ifstream tf_in(tf_file, std::ios::binary);
        if (!tf_in) {
            std::cerr << "Failed to open tf file: " << tf_file << std::endl;
            exit(1);
        }
        
        // Use smaller data type for loading to save memory
        std::vector<uint64_t> tf_values(TOTAL_13MERS);
        tf_in.read(reinterpret_cast<char*>(tf_values.data()), TOTAL_13MERS * sizeof(uint64_t));
        tf_in.close();
        
        // Allocate indices array with ARM64 alignment
        std::cout << "[ARM64] Allocating indices array..." << std::endl;
        
#ifdef __APPLE__
        void* aligned_mem = nullptr;
        if (posix_memalign(&aligned_mem, CACHE_LINE_SIZE, (TOTAL_13MERS + 1) * sizeof(uint64_t)) != 0) {
            throw std::bad_alloc();
        }
        indices = static_cast<uint64_t*>(aligned_mem);
#else
        indices = new uint64_t[TOTAL_13MERS + 1];
#endif
        
        if (indices == nullptr) {
            std::cerr << "Failed to allocate memory for indices" << std::endl;
            exit(1);
        }
        
        // Build cumulative indices efficiently
        indices[0] = 0;
        for (uint64_t i = 1; i < TOTAL_13MERS + 1; ++i) {
            indices[i] = indices[i-1] + tf_values[i-1];
            total_size += tf_values[i-1];
            max_tf = std::max(max_tf, tf_values[i-1]);
        }
        
        std::cout << "\t[ARM64] max_tf: " << max_tf << std::endl;
        std::cout << "\t[ARM64] total_size: " << total_size << std::endl;
        
        // Allocate position completion counters (32-bit to save memory)
        std::cout << "[ARM64] Allocating position completion counters..." << std::endl;
        
#ifdef __APPLE__
        if (posix_memalign(&aligned_mem, CACHE_LINE_SIZE, TOTAL_13MERS * sizeof(std::atomic<uint32_t>)) != 0) {
            throw std::bad_alloc();
        }
        ppositions = static_cast<std::atomic<uint32_t>*>(aligned_mem);
#else
        ppositions = new std::atomic<uint32_t>[TOTAL_13MERS];
#endif
        
        if (ppositions == nullptr) {
            std::cerr << "Failed to allocate memory for ppositions" << std::endl;
            exit(1);
        }
        
        // Zero-initialize efficiently
        memset(ppositions, 0, TOTAL_13MERS * sizeof(std::atomic<uint32_t>));
        
        // Allocate positions array
        std::cout << "[ARM64] Allocating positions array..." << std::endl;
        
#ifdef __APPLE__
        if (posix_memalign(&aligned_mem, CACHE_LINE_SIZE, total_size * sizeof(std::atomic<uint64_t>)) != 0) {
            throw std::bad_alloc();
        }
        positions = static_cast<std::atomic<uint64_t>*>(aligned_mem);
#else
        positions = new std::atomic<uint64_t>[total_size];
#endif
        
        if (positions == nullptr) {
            std::cerr << "Failed to allocate memory for positions" << std::endl;
            exit(1);
        }
        
        // Zero-initialize efficiently
        memset(positions, 0, total_size * sizeof(std::atomic<uint64_t>));
        
        std::cout << "[ARM64] AIndex13 structure initialized with optimizations" << std::endl;
    }
    
    ~ARM64AIndex13() {
#ifdef __APPLE__
        if (indices) free(indices);
        if (ppositions) free(ppositions);
        if (positions) free(positions);
#else
        delete[] indices;
        delete[] ppositions;
        delete[] positions;
#endif
    }
    
    /**
     * ARM64-optimized nucleotide validation
     */
    inline bool is_valid_kmer_fast(const char* kmer, int k) const {
        for (int i = 0; i < k; ++i) {
            if (!valid_nucleotide[static_cast<unsigned char>(kmer[i])]) {
                return false;
            }
        }
        return true;
    }
    
    /**
     * ARM64-optimized reverse complement using lookup table
     */
    std::string get_reverse_complement_fast(const char* kmer, int k) {
        std::string rc;
        rc.reserve(k);
        
        for (int i = k - 1; i >= 0; --i) {
            rc += static_cast<char>(complement_table[static_cast<unsigned char>(kmer[i])]);
        }
        
        return rc;
    }
    
    /**
     * ARM64-optimized index building with batch processing
     */
    void fill_index_from_reads(char* contents, uint64_t length, uint32_t num_threads, 
                              const HASHER& hasher) {
        std::cout << "[ARM64] Building 13-mer index from reads..." << std::endl;
        std::cout << "[ARM64] Length: " << length << ", Threads: " << num_threads << std::endl;
        
        // Optimize thread count for Apple Silicon
        if (num_threads > 8) {
            num_threads = 8; // Limit to P-cores
            std::cout << "[ARM64] Limited threads to " << num_threads << " for efficiency" << std::endl;
        }
        
        uint64_t batch_size = (length / num_threads) + 1;
        std::vector<std::thread> threads;
        
        for (uint32_t worker_id = 0; worker_id < num_threads; ++worker_id) {
            uint64_t start = worker_id * batch_size;
            uint64_t end = std::min((worker_id + 1) * batch_size, length);
            
            // Adjust start to avoid cutting k-mers
            if (start > 13) {
                start -= 12;
            }
            
            threads.emplace_back(&ARM64AIndex13::worker_thread, this, worker_id, start, end, 
                               contents, std::ref(hasher));
        }
        
        for (auto& t : threads) {
            t.join();
        }
        
        std::cout << "[ARM64] Index building completed" << std::endl;
    }
    
    /**
     * ARM64-optimized worker thread with batching
     */
    void worker_thread(int worker_id, uint64_t start, uint64_t end, char* contents, 
                      const HASHER& hasher) {
        emphf::stl_string_adaptor str_adapter;
        static std::mutex progress_mutex;
        
        const int k = 13;
        char ckmer[k + 1];
        ckmer[k] = '\0';
        
        uint64_t original_start = start;
        uint64_t original_end = end;
        
        // Skip to valid start position efficiently
        while (start < end - k + 1) {
            bool found = false;
            for (uint64_t i = start; i < start + k; ++i) {
                char c = contents[i];
                if (c == '\n' || c == '~' || c == '?' || c == 'N') {
                    start = i + 1;
                    found = true;
                    break;
                }
            }
            if (!found) {
                break;
            }
        }
        
        uint64_t total = (end > start + k - 1) ? (end - start - k + 1) : 0;
        uint64_t processed = 0;
        uint64_t valid_kmers = 0;
        uint64_t batch_count = 0;
        
        std::unique_lock<std::mutex> lock2(progress_mutex);
        std::cout << "[ARM64] Worker " << worker_id << " started, range: " << original_start 
                 << "-" << original_end << ", adjusted: " << start << "-" << end 
                 << ", total k-mers: " << total << std::endl;
        lock2.unlock();
        
        // Process in batches for better cache performance
        for (uint64_t i = start; i < end - k + 1; i += BATCH_SIZE) {
            uint64_t batch_end = std::min(i + BATCH_SIZE, end - k + 1);
            
            for (uint64_t j = i; j < batch_end; ++j) {
                processed++;
                
                // Fast nucleotide validation using lookup table
                bool skip = false;
                for (int l = 0; l < k; ++l) {
                    char c = contents[j + l];
                    if (c == '\n' || c == '~' || c == '?' || !valid_nucleotide[static_cast<unsigned char>(c)]) {
                        skip = true;
                        break;
                    }
                }
                
                if (skip) {
                    continue;
                }
                
                // Extract k-mer efficiently
                std::memcpy(ckmer, &contents[j], k);
                
                // Process k-mer with minimal string operations
                bool stored = false;
                
                // Try forward direction
                try {
                    std::string kmer(ckmer, k);
                    uint64_t hash_id = hasher.lookup(kmer, str_adapter);
                    if (hash_id < TOTAL_13MERS) {
                        uint32_t pos_idx = ppositions[hash_id].fetch_add(1, std::memory_order_relaxed);
                        uint64_t array_idx = indices[hash_id] + pos_idx;
                        
                        if (array_idx < total_size && pos_idx < (indices[hash_id + 1] - indices[hash_id])) {
                            positions[array_idx].store(j + 1, std::memory_order_relaxed); // 1-based position
                            valid_kmers++;
                            stored = true;
                        }
                    }
                } catch (...) {
                    // Forward lookup failed
                }
                
                // Try reverse complement if forward failed
                if (!stored) {
                    try {
                        std::string rev_kmer = get_reverse_complement_fast(ckmer, k);
                        uint64_t hash_id = hasher.lookup(rev_kmer, str_adapter);
                        if (hash_id < TOTAL_13MERS) {
                            uint32_t pos_idx = ppositions[hash_id].fetch_add(1, std::memory_order_relaxed);
                            uint64_t array_idx = indices[hash_id] + pos_idx;
                            
                            if (array_idx < total_size && pos_idx < (indices[hash_id + 1] - indices[hash_id])) {
                                positions[array_idx].store(j + 1, std::memory_order_relaxed);
                                valid_kmers++;
                            }
                        }
                    } catch (...) {
                        // Both failed - skip
                    }
                }
            }
            
            // Progress reporting per batch
            batch_count++;
            if (batch_count % 100 == 0) {
                int progress = total > 0 ? int((100.0 * processed) / total) : 100;
                std::unique_lock<std::mutex> lock(progress_mutex);
                std::cout << "[ARM64] Worker " << worker_id << " progress: " << progress 
                         << "%, processed: " << processed << "/" << total 
                         << ", valid k-mers: " << valid_kmers << std::endl;
            }
        }
        
        std::unique_lock<std::mutex> final_lock(progress_mutex);
        std::cout << "[ARM64] Worker " << worker_id << " finished. Valid k-mers: " << valid_kmers << std::endl;
    }
    
    /**
     * ARM64-optimized position retrieval
     */
    void get_positions(const std::string& kmer, const HASHER& hasher, std::vector<uint64_t>& result) {
        result.clear();
        result.reserve(64); // Pre-allocate for common case
        
        emphf::stl_string_adaptor str_adapter;
        
        // Try forward direction
        try {
            uint64_t hash_id = hasher.lookup(kmer, str_adapter);
            if (hash_id < TOTAL_13MERS) {
                uint64_t start_idx = indices[hash_id];
                uint64_t end_idx = indices[hash_id + 1];
                
                for (uint64_t i = start_idx; i < end_idx; ++i) {
                    uint64_t pos = positions[i].load(std::memory_order_relaxed);
                    if (pos > 0) { // Valid position
                        result.push_back(pos);
                    }
                }
                return;
            }
        } catch (...) {
            // Try reverse complement
            try {
                std::string rev_kmer = get_reverse_complement_fast(kmer.c_str(), kmer.length());
                uint64_t hash_id = hasher.lookup(rev_kmer, str_adapter);
                if (hash_id < TOTAL_13MERS) {
                    uint64_t start_idx = indices[hash_id];
                    uint64_t end_idx = indices[hash_id + 1];
                    
                    for (uint64_t i = start_idx; i < end_idx; ++i) {
                        uint64_t pos = positions[i].load(std::memory_order_relaxed);
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
    
    /**
     * ARM64-optimized save with efficient I/O
     */
    void save(const std::string& pos_bin_file, const std::string& index_bin_file, 
              const std::string& indices_bin_file, const std::vector<uint64_t>& start_positions) {
        
        std::cout << "[ARM64] Saving pos.bin array..." << std::endl;
        std::ofstream pos_out(pos_bin_file, std::ios::binary);
        pos_out.write(reinterpret_cast<const char*>(start_positions.data()), 
                     start_positions.size() * sizeof(uint64_t));
        pos_out.close();
        
        std::cout << "[ARM64] Saving index.bin array..." << std::endl;
        std::ofstream index_out(index_bin_file, std::ios::binary);
        std::cout << "[ARM64] Positions array size: " << sizeof(uint64_t) * total_size << std::endl;
        
        // Convert atomic to regular array for saving (efficient batch conversion)
        const size_t batch_size = 1024 * 1024; // 1M elements per batch
        std::vector<uint64_t> batch_buffer(batch_size);
        
        for (uint64_t i = 0; i < total_size; i += batch_size) {
            size_t current_batch_size = std::min(static_cast<size_t>(batch_size), static_cast<size_t>(total_size - i));
            
            for (size_t j = 0; j < current_batch_size; ++j) {
                batch_buffer[j] = positions[i + j].load(std::memory_order_relaxed);
            }
            
            index_out.write(reinterpret_cast<const char*>(batch_buffer.data()), 
                           sizeof(uint64_t) * current_batch_size);
        }
        
        index_out.close();
        
        std::cout << "[ARM64] Saving indices.bin array..." << std::endl;
        std::ofstream indices_out(indices_bin_file, std::ios::binary);
        indices_out.write(reinterpret_cast<const char*>(indices), 
                         sizeof(uint64_t) * (TOTAL_13MERS + 1));
        indices_out.close();
        
        std::cout << "[ARM64] All files saved successfully" << std::endl;
    }
};

// Static lookup tables for ARM64 optimization
const bool ARM64AIndex13::valid_nucleotide[256] = {
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, true,  false, true,  false, false, false, true,  false, false, false, false, false, false, false, false, // A, C, G
    false, false, false, false, true,  false, false, false, false, false, false, false, false, false, false, false, // T
    false, true,  false, true,  false, false, false, true,  false, false, false, false, false, false, false, false, // a, c, g
    false, false, false, false, true,  false, false, false, false, false, false, false, false, false, false, false, // t
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false
};

const unsigned char ARM64AIndex13::complement_table[256] = {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15,
    16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
    32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
    48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
    64,  'T', 'B', 'G', 'D', 'E', 'F', 'C', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
    'P', 'Q', 'R', 'S', 'A', 'U', 'V', 'W', 'X', 'Y', 'Z', 91,  92,  93,  94,  95,
    96,  'T', 'B', 'G', 'D', 'E', 'F', 'C', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
    'P', 'Q', 'R', 'S', 'A', 'U', 'V', 'W', 'X', 'Y', 'Z', 123, 124, 125, 126, 127,
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
    144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
    160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
    176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
    208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
    224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
    240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};

int main(int argc, char** argv) {
    if (argc < 6) {
        std::cerr << "ARM64-Optimized Compute AIndex for 13-mers with perfect hash" << std::endl;
        std::cerr << "Optimized for Apple Silicon (M1/M2/M3)" << std::endl;
        std::cerr << std::endl;
        std::cerr << "Usage: " << argv[0]
                 << " <reads_file> <hash_file> <tf_file> <output_prefix> <num_threads> [pos_bin] [index_bin] [indices_bin]" << std::endl;
        std::cerr << std::endl;
        std::cerr << "Arguments:" << std::endl;
        std::cerr << "  reads_file:      Input reads file (one sequence per line)" << std::endl;
        std::cerr << "  hash_file:       13-mer perfect hash file (.hash)" << std::endl;
        std::cerr << "  tf_file:         13-mer frequencies file (.tf.bin)" << std::endl;
        std::cerr << "  output_prefix:   Prefix for output files" << std::endl;
        std::cerr << "  num_threads:     Number of threads (auto-optimized for Apple Silicon)" << std::endl;
        std::cerr << "  pos_bin:         Optional: Output pos.bin filename" << std::endl;
        std::cerr << "  index_bin:       Optional: Output index.bin filename" << std::endl;
        std::cerr << "  indices_bin:     Optional: Output indices.bin filename" << std::endl;
        std::cerr << std::endl;
        std::cerr << "ARM64 Features:" << std::endl;
        std::cerr << "  - Cache-aligned data structures for better performance" << std::endl;
        std::cerr << "  - 32-bit position counters to reduce memory usage" << std::endl;
        std::cerr << "  - Batch processing for improved cache efficiency" << std::endl;
        std::cerr << "  - Optimized lookup tables for nucleotide validation" << std::endl;
        std::cerr << "  - Thread count auto-optimization for Apple Silicon P-cores" << std::endl;
        return 1;
    }

    std::string reads_file = argv[1];
    std::string hash_file = argv[2];
    std::string tf_file = argv[3];
    std::string output_prefix = argv[4];
    uint32_t num_threads = std::atoi(argv[5]);

    // Optional output filenames
    std::string pos_bin_file = (argc > 6) ? argv[6] : output_prefix + ".pos.bin";
    std::string index_bin_file = (argc > 7) ? argv[7] : output_prefix + ".index.bin";
    std::string indices_bin_file = (argc > 8) ? argv[8] : output_prefix + ".indices.bin";

    std::cout << "=== ARM64-Optimized AIndex13 Computation ===" << std::endl;
    std::cout << "Platform: Apple Silicon optimized" << std::endl;
    std::cout << "Reads file: " << reads_file << std::endl;
    std::cout << "Hash file: " << hash_file << std::endl;
    std::cout << "TF file: " << tf_file << std::endl;
    std::cout << "Output prefix: " << output_prefix << std::endl;
    std::cout << "Threads: " << num_threads << " (will be optimized for ARM64)" << std::endl;
    std::cout << std::endl;

    std::cout << "[ARM64] Loading 13-mer hash..." << std::endl;
    
    // Load hash
    HASHER hasher;
    std::ifstream hash_in(hash_file, std::ios::binary);
    if (!hash_in) {
        std::cerr << "Failed to open hash file: " << hash_file << std::endl;
        return 1;
    }
    hasher.load(hash_in);
    hash_in.close();
    
    std::cout << "[ARM64] Hash loaded successfully" << std::endl;

    std::cout << "[ARM64] Loading reads and building index..." << std::endl;
    
    // Load reads file efficiently
    std::ifstream reads_in(reads_file);
    if (!reads_in) {
        std::cerr << "Failed to open reads file: " << reads_file << std::endl;
        return 1;
    }

    reads_in.seekg(0, std::ios::end);
    uint64_t length = reads_in.tellg();
    
    char* contents = nullptr;
#ifdef __APPLE__
    // Use aligned allocation for better ARM64 performance
    if (posix_memalign(reinterpret_cast<void**>(&contents), 4096, length + 1) != 0) {
        std::cerr << "Failed to allocate aligned memory for reads: " << length + 1 << std::endl;
        return 1;
    }
#else
    contents = new char[length + 1];
    if (contents == nullptr) {
        std::cerr << "Failed to allocate memory for reads: " << length + 1 << std::endl;
        return 1;
    }
#endif

    reads_in.seekg(0, std::ios::beg);
    reads_in.read(contents, length);
    reads_in.close();
    contents[length] = 0;
    
    // Build read start positions efficiently
    std::vector<uint64_t> start_positions;
    start_positions.reserve(length / 50); // Reasonable estimate
    std::unordered_map<uint64_t, uint32_t> start2rid;
    
    uint32_t rid = 0;
    uint64_t pos = 0;
    start_positions.push_back(pos);
    start2rid[pos] = rid;
    
    uint64_t nreads = 0;
    for (uint64_t i = 0; i < length; i++) {
        if (contents[i] == '\n') {
            start_positions.push_back(i + 1);
            rid++;
            start2rid[i + 1] = rid;
            nreads++;
            
            if (rid % 1000000 == 0) {
                std::cout << "[ARM64] Loaded reads: " << rid << std::endl;
            }
        }
    }
    
    start_positions.push_back(length);
    std::cout << "[ARM64] Loaded " << nreads << " reads, " << length << " characters" << std::endl;

    // Initialize ARM64AIndex13 and build index
    std::cout << "[ARM64] Initializing AIndex13..." << std::endl;
    ARM64AIndex13 aindex(tf_file);
    
    std::cout << "[ARM64] Building index from reads..." << std::endl;
    aindex.fill_index_from_reads(contents, length, num_threads, hasher);
    
    std::cout << "[ARM64] Saving index..." << std::endl;
    aindex.save(pos_bin_file, index_bin_file, indices_bin_file, start_positions);

#ifdef __APPLE__
    free(contents);
#else
    delete[] contents;
#endif
    
    std::cout << "[ARM64] AIndex13 computation completed successfully!" << std::endl;
    return 0;
}
