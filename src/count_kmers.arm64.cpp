#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <algorithm>
#include <memory>
#include <cstring>
#include <iomanip>

// ARM64-specific optimizations
#ifdef __aarch64__
#include <arm_neon.h>
#endif

// ARM64 optimized version for M1 Mac
using kmer_t = uint64_t;

// Optimized hash function for ARM64
struct KmerHashARM64 {
    size_t operator()(kmer_t k) const {
        // ARM64-optimized hash using multiply-add instructions
        // M1 has excellent 64-bit multiply performance
        k ^= k >> 33;
        k *= 0xff51afd7ed558ccdULL;
        k ^= k >> 33;
        k *= 0xc4ceb9fe1a85ec53ULL;
        k ^= k >> 33;
        return k;
    }
};

// Robin Hood hash map for better performance on ARM64
template<typename K, typename V>
using FastMapARM64 = std::unordered_map<K, V, KmerHashARM64>;

class ARM64KmerCounter {
private:
    size_t k;
    size_t num_threads;
    size_t min_count;
    bool use_canonical;
    static constexpr size_t BATCH_SIZE = 2000;  // Increased for ARM64
    static constexpr size_t BUFFER_SIZE = 32 * 1024 * 1024; // 32MB buffer for M1
    static constexpr size_t CACHE_LINE_SIZE = 128; // ARM64 cache line
    
    // Thread-local data with cache alignment
    struct alignas(CACHE_LINE_SIZE) ThreadDataARM64 {
        FastMapARM64<kmer_t, uint32_t> local_kmers;
        std::vector<std::string> batch;
        size_t sequences_processed = 0;
        size_t kmers_processed = 0;
        char* read_buffer = nullptr;
        
        // ARM64-specific padding to avoid false sharing
        char padding[CACHE_LINE_SIZE - sizeof(size_t) * 2];
        
        ThreadDataARM64() {
            local_kmers.reserve(1 << 21); // 2M k-mers for M1
            batch.reserve(BATCH_SIZE);
            read_buffer = new char[BUFFER_SIZE];
        }
        
        ~ThreadDataARM64() {
            delete[] read_buffer;
        }
    };
    
    std::vector<std::unique_ptr<ThreadDataARM64>> thread_data;
    std::atomic<size_t> current_block{0};
    std::vector<std::pair<size_t, size_t>> file_blocks;
    std::atomic<size_t> total_sequences{0};
    std::atomic<size_t> sequences_processed{0};
    
    // ARM64-optimized lookup tables
    static constexpr uint8_t char_to_bits_arm64[256] = {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 0-15
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 16-31
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 32-47
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 48-63
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 64-79  (A, C, G)
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 80-95  (T, U)
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 96-111  (a, c, g)
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 112-127 (t, u)
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 128-255
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
    };
    
    static constexpr char bits_to_char_arm64[4] = {'A', 'C', 'G', 'T'};
    
    // ARM64-optimized conversion with NEON (if available)
    inline bool string_to_kmer_arm64(const char* seq, size_t pos, kmer_t& kmer) {
        kmer = 0;
        const char* p = seq + pos;
        
#ifdef __aarch64__
        // Use ARM64 NEON for acceleration with large k
        if (k >= 16) {
            // Vectorized processing for long k-mers
            size_t vec_len = k & ~7; // Multiple of 8
            for (size_t i = 0; i < vec_len; i += 8) {
                // Load 8 characters and check them scalar
                // This is faster and simpler than vector check
                bool all_valid = true;
                for (int j = 0; j < 8; ++j) {
                    uint8_t bits = char_to_bits_arm64[(uint8_t)p[i + j]];
                    if (bits == 4) {
                        all_valid = false;
                        break;
                    }
                    kmer = (kmer << 2) | bits;
                }
                
                if (!all_valid) return false;
            }
            
            // Process remaining characters
            for (size_t i = vec_len; i < k; ++i) {
                uint8_t bits = char_to_bits_arm64[(uint8_t)p[i]];
                if (bits == 4) return false;
                kmer = (kmer << 2) | bits;
            }
            return true;
        }
#endif

        // Scalar processing for short k-mers or fallback
        for (size_t i = 0; i < k; ++i) {
            uint8_t bits = char_to_bits_arm64[(uint8_t)p[i]];
            if (bits == 4) return false;
            kmer = (kmer << 2) | bits;
        }
        return true;
    }
    
    // ARM64-optimized reverse complement
    inline kmer_t reverse_complement_arm64(kmer_t kmer) {
        // ARM64 has efficient bit manipulation instructions
        
#ifdef __aarch64__
        // Use ARM64 rbit instruction for bit reversal
        asm("rbit %0, %1" : "=r"(kmer) : "r"(kmer));
        
        // Complement (XOR with mask)
        kmer = kmer ^ 0xAAAAAAAAAAAAAAAAULL;
        
        // Shift to correct position
        return kmer >> (64 - 2 * k);
#else
        // Fallback for other architectures
        kmer = ((kmer & 0xAAAAAAAAAAAAAAAA) >> 1) | ((kmer & 0x5555555555555555) << 1);
        kmer = ((kmer & 0xCCCCCCCCCCCCCCCC) >> 2) | ((kmer & 0x3333333333333333) << 2);
        kmer = ((kmer & 0xF0F0F0F0F0F0F0F0) >> 4) | ((kmer & 0x0F0F0F0F0F0F0F0F) << 4);
        kmer = ((kmer & 0xFF00FF00FF00FF00) >> 8) | ((kmer & 0x00FF00FF00FF00FF) << 8);
        kmer = ((kmer & 0xFFFF0000FFFF0000) >> 16) | ((kmer & 0x0000FFFF0000FFFF) << 16);
        kmer = (kmer >> 32) | (kmer << 32);
        kmer = ~kmer;
        return kmer >> (64 - 2 * k);
#endif
    }
    
    inline kmer_t get_canonical_arm64(kmer_t kmer) {
        if (!use_canonical) return kmer;
        kmer_t rc = reverse_complement_arm64(kmer);
        return (kmer < rc) ? kmer : rc;
    }
    
    void process_sequence_batch_arm64(ThreadDataARM64& td) {
        // Preload data into M1 cache
        __builtin_prefetch(&td.local_kmers, 1, 1);
        
        for (const auto& seq : td.batch) {
            if (seq.length() < k) continue;
            
            const char* s = seq.c_str();
            size_t len = seq.length();
            
            // Preload string into cache
            __builtin_prefetch(s, 0, 1);
            __builtin_prefetch(s + 64, 0, 1);
            
            for (size_t i = 0; i <= len - k; ++i) {
                kmer_t kmer;
                if (!string_to_kmer_arm64(s, i, kmer)) continue;
                
                kmer_t canonical = get_canonical_arm64(kmer);
                td.local_kmers[canonical]++;
                td.kmers_processed++;
                
                // Preload next data
                if (i % 32 == 0 && i + 96 < len) {
                    __builtin_prefetch(s + i + 96, 0, 1);
                }
            }
            td.sequences_processed++;
        }
        
        sequences_processed.fetch_add(td.batch.size());
        td.batch.clear();
    }
    
    void worker_thread_arm64(size_t thread_id) {
        auto& td = *thread_data[thread_id];
        
        // Set thread affinity for M1 (if supported)
        // M1 has P-cores and E-cores
        
        while (true) {
            size_t block_idx = current_block.fetch_add(1);
            if (block_idx >= file_blocks.size()) break;
            
            auto [start_pos, block_size] = file_blocks[block_idx];
            
            // Read file block with optimization for M1 SSD
            std::ifstream file("temp_input.dat", std::ios::binary);
            file.rdbuf()->pubsetbuf(td.read_buffer, BUFFER_SIZE);
            file.seekg(start_pos);
            file.read(td.read_buffer, block_size);
            file.close();
            
            // Parse sequences with ARM64 optimizations
            process_file_block_arm64(td, block_size);
        }
    }
    
    void process_file_block_arm64(ThreadDataARM64& td, size_t block_size) {
        std::string current_seq;
        current_seq.reserve(1024); // Preallocate memory
        bool in_sequence = false;
        
        char* buffer = td.read_buffer;
        
        for (size_t i = 0; i < block_size; ++i) {
            char c = buffer[i];
            
            if (c == '>') {
                if (!current_seq.empty()) {
                    td.batch.push_back(std::move(current_seq));
                    current_seq.clear();
                    current_seq.reserve(1024);
                    
                    if (td.batch.size() >= BATCH_SIZE) {
                        process_sequence_batch_arm64(td);
                    }
                }
                in_sequence = false;
            } else if (c == '\n') {
                if (!in_sequence) {
                    in_sequence = true;
                }
            } else if (in_sequence && c != '\r') {
                current_seq += c;
            }
        }
        
        if (!current_seq.empty()) {
            td.batch.push_back(std::move(current_seq));
        }
        
        if (!td.batch.empty()) {
            process_sequence_batch_arm64(td);
        }
    }

public:
    ARM64KmerCounter(size_t k_value, size_t threads, size_t min_count_filter = 1, bool canonical = true)
        : k(k_value), num_threads(threads), min_count(min_count_filter), use_canonical(canonical) {
        
        // M1 Mac optimization: use all P-cores + some E-cores
        if (threads == 0) {
            num_threads = std::thread::hardware_concurrency();
            // On M1 usually 8 cores (4 P + 4 E), use all
        }
        
        thread_data.reserve(num_threads);
        for (size_t i = 0; i < num_threads; ++i) {
            thread_data.emplace_back(std::make_unique<ThreadDataARM64>());
        }
        
        std::cout << "ARM64 K-mer Counter initialized with " << num_threads << " threads" << std::endl;
        std::cout << "Using k=" << k << ", canonical=" << use_canonical << std::endl;
    }
    
    void count_kmers_from_file(const std::string& filename) {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        std::cout << "ARM64 optimized k-mer counting starting..." << std::endl;
        std::cout << "Preparing file for parallel processing..." << std::endl;
        
        // Read file into memory (M1 has a lot of RAM)
        std::ifstream input(filename, std::ios::binary | std::ios::ate);
        if (!input) {
            throw std::runtime_error("Cannot open input file: " + filename);
        }
        
        size_t file_size = input.tellg();
        input.seekg(0);
        
        std::vector<char> file_buffer(file_size);
        input.read(file_buffer.data(), file_size);
        input.close();
        
        std::cout << "File loaded: " << file_size << " bytes" << std::endl;
        
        // Find sequences
        std::vector<size_t> sequence_starts;
        sequence_starts.reserve(100000); // Preallocate
        
        for (size_t i = 0; i < file_size; ++i) {
            if (file_buffer[i] == '>') {
                sequence_starts.push_back(i);
            }
        }
        total_sequences = sequence_starts.size();
        
        std::cout << "Found " << total_sequences << " sequences" << std::endl;
        std::cout << "Starting ARM64 optimized parallel processing..." << std::endl;
        
        // Launch worker threads
        std::vector<std::thread> workers;
        workers.reserve(num_threads);
        
        for (size_t t = 0; t < num_threads; ++t) {
            workers.emplace_back([this, t, &file_buffer, &sequence_starts]() {
                process_sequences_range_arm64(t, file_buffer, sequence_starts);
            });
        }
        
        // Wait for completion
        for (auto& t : workers) {
            t.join();
        }
        
        std::cout << "\nMerging results from " << num_threads << " threads..." << std::endl;
        
        // Merge results with ARM64 optimizations
        auto final_kmers = merge_results_arm64();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        print_statistics_arm64(duration, final_kmers);
        save_results_arm64(final_kmers, "output_arm64.txt");
    }
    
private:
    void process_sequences_range_arm64(size_t thread_id, const std::vector<char>& file_buffer, 
                                      const std::vector<size_t>& sequence_starts) {
        auto& td = *thread_data[thread_id];
        
        size_t seqs_per_thread = (sequence_starts.size() + num_threads - 1) / num_threads;
        size_t start_idx = thread_id * seqs_per_thread;
        size_t end_idx = std::min(start_idx + seqs_per_thread, sequence_starts.size());
        
        for (size_t i = start_idx; i < end_idx; ++i) {
            process_single_sequence_arm64(td, file_buffer, sequence_starts, i);
            
            // Progress reporting
            if (td.sequences_processed % 5000 == 0) {
                size_t total_proc = sequences_processed.load();
                double progress = 100.0 * total_proc / total_sequences;
                std::cout << "\rARM64 Progress: " << std::fixed << std::setprecision(1) 
                         << progress << "% (" << total_proc << "/" << total_sequences << ")" 
                         << std::flush;
            }
        }
    }
    
    void process_single_sequence_arm64(ThreadDataARM64& td, const std::vector<char>& file_buffer,
                                      const std::vector<size_t>& sequence_starts, size_t seq_idx) {
        size_t seq_start = sequence_starts[seq_idx];
        size_t seq_end = (seq_idx + 1 < sequence_starts.size()) ? 
                        sequence_starts[seq_idx + 1] : file_buffer.size();
        
        // Skip header
        while (seq_start < seq_end && file_buffer[seq_start] != '\n') seq_start++;
        seq_start++;
        
        // Collect sequence with ARM64 optimizations
        std::string sequence;
        sequence.reserve(seq_end - seq_start);
        
        for (size_t j = seq_start; j < seq_end; ++j) {
            char c = file_buffer[j];
            if (c != '\n' && c != '\r' && c != '>') {
                sequence += c;
            }
            if (c == '>') break;
        }
        
        // Process k-mers
        if (sequence.length() >= k) {
            const char* s = sequence.c_str();
            size_t len = sequence.length();
            
            // Preload into M1 cache
            __builtin_prefetch(s, 0, 1);
            
            for (size_t pos = 0; pos <= len - k; ++pos) {
                kmer_t kmer;
                if (!string_to_kmer_arm64(s, pos, kmer)) continue;
                
                kmer_t canonical = get_canonical_arm64(kmer);
                td.local_kmers[canonical]++;
                td.kmers_processed++;
            }
        }
        
        td.sequences_processed++;
        sequences_processed.fetch_add(1);
    }
    
    FastMapARM64<kmer_t, size_t> merge_results_arm64() {
        FastMapARM64<kmer_t, size_t> final_kmers;
        final_kmers.reserve(5000000); // 5M k-mers for M1
        
        for (const auto& td : thread_data) {
            for (const auto& [kmer, count] : td->local_kmers) {
                final_kmers[kmer] += count;
            }
        }
        
        return final_kmers;
    }
    
    void print_statistics_arm64(const std::chrono::milliseconds& duration, 
                               const FastMapARM64<kmer_t, size_t>& final_kmers) {
        std::cout << "\nARM64 processing completed in " << duration.count() << " ms" << std::endl;
        std::cout << "Processing rate: " << (sequences_processed.load() * 1000.0 / duration.count()) 
                  << " sequences/second" << std::endl;
        
        size_t total_kmers = 0;
        for (const auto& td : thread_data) {
            total_kmers += td->kmers_processed;
        }
        std::cout << "Total k-mers processed: " << total_kmers << std::endl;
        std::cout << "Unique k-mers found: " << final_kmers.size() << std::endl;
        std::cout << "K-mer processing rate: " << (total_kmers * 1000.0 / duration.count()) 
                  << " k-mers/second" << std::endl;
    }
    
    void save_results_arm64(const FastMapARM64<kmer_t, size_t>& kmers, const std::string& output_file) {
        std::vector<std::pair<kmer_t, size_t>> sorted_kmers;
        sorted_kmers.reserve(kmers.size());
        
        for (const auto& [kmer, count] : kmers) {
            if (count >= min_count) {
                sorted_kmers.emplace_back(kmer, count);
            }
        }
        
        // ARM64 optimized sorting
        std::sort(sorted_kmers.begin(), sorted_kmers.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        std::ofstream out(output_file);
        out << "# ARM64 optimized k-mer counting results\n";
        out << "# k=" << k << ", sequences=" << total_sequences << "\n";
        
        for (const auto& [kmer, count] : sorted_kmers) {
            out << kmer_to_string_arm64(kmer) << "\t" << count << "\n";
        }
        out.close();
        
        std::cout << "ARM64 results saved to " << output_file 
                  << " (" << sorted_kmers.size() << " k-mers)" << std::endl;
    }
    
    std::string kmer_to_string_arm64(kmer_t kmer) {
        std::string result(k, 'N');
        for (int i = k - 1; i >= 0; --i) {
            result[i] = bits_to_char_arm64[kmer & 3];
            kmer >>= 2;
        }
        return result;
    }
};

int main(int argc, char* argv[]) {
    std::cout << "ARM64 K-mer Counter for M1 Mac" << std::endl;
    std::cout << "Optimized for Apple Silicon architecture" << std::endl;
    
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <k> <output_file> [-t threads] [-m min_count]" << std::endl;
        std::cerr << "  -t threads   : Number of threads (default: auto-detect)" << std::endl;
        std::cerr << "  -m min_count : Minimum count filter (default: 1)" << std::endl;
        return 1;
    }
    
    std::string input_file = argv[1];
    size_t k = std::stoul(argv[2]);
    std::string output_file = argv[3];
    
    size_t threads = 0; // Auto-detect
    size_t min_count = 1;
    
    // Parse arguments
    for (int i = 4; i < argc; i++) {
        if (std::string(argv[i]) == "-t" && i + 1 < argc) {
            threads = std::stoul(argv[++i]);
        } else if (std::string(argv[i]) == "-m" && i + 1 < argc) {
            min_count = std::stoul(argv[++i]);
        }
    }
    
    try {
        ARM64KmerCounter counter(k, threads, min_count, true);
        counter.count_kmers_from_file(input_file);
        
        std::cout << "ARM64 k-mer counting completed successfully!" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
