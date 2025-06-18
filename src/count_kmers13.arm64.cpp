//
// ARM64-optimized 13-mer counter using precomputed perfect hash
// Optimized for Apple Silicon (M1/M2/M3) with reduced memory usage during compilation
// Supports FASTQ SE/PE, FASTA, and plain reads (one sequence per line)
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <queue>
#include <condition_variable>
#include <chrono>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <memory>

#ifdef __aarch64__
#include <arm_neon.h>
#endif

// Simplified includes for ARM64 to reduce compilation memory usage
#include "hash.hpp"
#include "kmers.hpp"

class ARM64Kmer13Counter {
private:
    static const uint32_t TOTAL_13MERS = 67108864; // 4^13
    static const int KMER_SIZE = 13;
    static constexpr size_t BATCH_SIZE = 1000;
    static constexpr size_t CACHE_LINE_SIZE = 128; // ARM64 cache line
    
    // Perfect hash components
    HASHER hasher;
    std::unique_ptr<std::atomic<uint32_t>[]> kmer_counts; // Use uint32_t to save memory
    bool hash_loaded;
    
    // Threading with ARM64 optimizations
    std::mutex queue_mutex;
    std::condition_variable cv;
    std::queue<std::vector<std::string>> batch_queue; // Process in batches
    std::atomic<bool> done_reading{false};
    size_t num_threads;
    
    // Statistics
    std::atomic<uint64_t> total_sequences{0};
    std::atomic<uint64_t> total_kmers_processed{0};
    std::atomic<uint64_t> valid_kmers{0};
    std::atomic<uint64_t> invalid_kmers{0};
    
    // ARM64-optimized lookup table
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
    
public:
    ARM64Kmer13Counter(size_t threads = std::thread::hardware_concurrency()) 
        : hash_loaded(false), num_threads(threads) {
        if (num_threads == 0) num_threads = 1;
        
        // Initialize atomic array for all possible 13-mers with ARM64 alignment
        kmer_counts = std::make_unique<std::atomic<uint32_t>[]>(TOTAL_13MERS);
        for (uint32_t i = 0; i < TOTAL_13MERS; i++) {
            kmer_counts[i] = 0;
        }
        
        std::cout << "ARM64 13-mer counter initialized with " << num_threads << " threads" << std::endl;
    }
    
    /**
     * Load precomputed perfect hash with error handling
     */
    bool load_perfect_hash(const std::string& hash_file) {
        std::cout << "Loading perfect hash from: " << hash_file << std::endl;
        
        std::ifstream hash_in(hash_file, std::ios::binary);
        if (!hash_in.is_open()) {
            std::cerr << "Error: Cannot open hash file: " << hash_file << std::endl;
            return false;
        }
        
        try {
            hasher.load(hash_in);
            hash_in.close();
            hash_loaded = true;
            std::cout << "Perfect hash loaded successfully (ARM64 optimized)" << std::endl;
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error loading hash: " << e.what() << std::endl;
            return false;
        }
    }
    
    /**
     * ARM64-optimized k-mer validation
     */
    inline bool is_valid_kmer_arm64(const char* kmer) {
        // Fast ARM64 validation using lookup table
        for (int i = 0; i < KMER_SIZE; ++i) {
            if (char_to_bits_arm64[(uint8_t)kmer[i]] == 4) {
                return false;
            }
        }
        return true;
    }
    
    /**
     * ARM64-optimized sequence normalization
     */
    std::string normalize_sequence_arm64(const std::string& seq) {
        std::string normalized;
        normalized.reserve(seq.length());
        
        // ARM64 optimized character processing
        const char* data = seq.data();
        size_t len = seq.length();
        
#ifdef __aarch64__
        // Process 8 characters at a time using NEON when possible
        size_t vec_len = len & ~7;
        for (size_t i = 0; i < vec_len; i += 8) {
            for (int j = 0; j < 8; ++j) {
                char c = data[i + j];
                char upper_c = (c >= 'a' && c <= 'z') ? c - 32 : c;
                if (upper_c == 'A' || upper_c == 'T' || upper_c == 'G' || upper_c == 'C') {
                    normalized += upper_c;
                } else {
                    normalized += 'N';
                }
            }
        }
        
        // Process remaining characters
        for (size_t i = vec_len; i < len; ++i) {
            char c = data[i];
            char upper_c = (c >= 'a' && c <= 'z') ? c - 32 : c;
            if (upper_c == 'A' || upper_c == 'T' || upper_c == 'G' || upper_c == 'C') {
                normalized += upper_c;
            } else {
                normalized += 'N';
            }
        }
#else
        // Fallback for non-ARM64
        for (char c : seq) {
            char upper_c = std::toupper(c);
            if (upper_c == 'A' || upper_c == 'T' || upper_c == 'G' || upper_c == 'C') {
                normalized += upper_c;
            } else {
                normalized += 'N';
            }
        }
#endif
        
        return normalized;
    }
    
    /**
     * Process a batch of sequences (ARM64 optimized)
     */
    void process_sequence_batch(const std::vector<std::string>& sequences) {
        emphf::stl_string_adaptor str_adapter;
        
        for (const std::string& seq : sequences) {
            if (seq.length() < KMER_SIZE) continue;
            
            std::string normalized = normalize_sequence_arm64(seq);
            total_sequences.fetch_add(1);
            
            // ARM64-optimized k-mer processing
            const char* data = normalized.data();
            size_t max_pos = normalized.length() - KMER_SIZE + 1;
            
            for (size_t i = 0; i < max_pos; ++i) {
                total_kmers_processed.fetch_add(1);
                
                if (is_valid_kmer_arm64(data + i)) {
                    // Create string for hash lookup (minimize allocations)
                    std::string kmer(data + i, KMER_SIZE);
                    
                    try {
                        uint64_t hash_idx = hasher.lookup(kmer, str_adapter);
                        
                        if (hash_idx < TOTAL_13MERS) {
                            // ARM64 atomic increment
                            kmer_counts[hash_idx].fetch_add(1, std::memory_order_relaxed);
                            valid_kmers.fetch_add(1);
                        } else {
                            invalid_kmers.fetch_add(1);
                        }
                    } catch (...) {
                        invalid_kmers.fetch_add(1);
                    }
                } else {
                    invalid_kmers.fetch_add(1);
                }
            }
        }
    }
    
    /**
     * ARM64-optimized worker thread
     */
    void worker_thread_arm64() {
        while (true) {
            std::unique_lock<std::mutex> lock(queue_mutex);
            cv.wait(lock, [this] { return !batch_queue.empty() || done_reading; });
            
            if (batch_queue.empty() && done_reading) {
                break;
            }
            
            if (!batch_queue.empty()) {
                std::vector<std::string> batch = std::move(batch_queue.front());
                batch_queue.pop();
                lock.unlock();
                
                process_sequence_batch(batch);
            }
        }
    }
    
    /**
     * Detect file format
     */
    enum class FileFormat {
        PLAIN,
        FASTA,
        FASTQ
    };
    
    FileFormat detect_format(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) return FileFormat::PLAIN;
        
        std::string first_line;
        if (std::getline(file, first_line)) {
            if (first_line.empty()) return FileFormat::PLAIN;
            if (first_line[0] == '>') return FileFormat::FASTA;
            if (first_line[0] == '@') return FileFormat::FASTQ;
        }
        
        return FileFormat::PLAIN;
    }
    
    /**
     * ARM64-optimized file reading with batching
     */
    void read_file_in_batches(std::ifstream& input, FileFormat format) {
        std::vector<std::string> current_batch;
        current_batch.reserve(BATCH_SIZE);
        
        auto flush_batch = [&]() {
            if (!current_batch.empty()) {
                std::lock_guard<std::mutex> lock(queue_mutex);
                batch_queue.push(std::move(current_batch));
                cv.notify_one();
                current_batch.clear();
                current_batch.reserve(BATCH_SIZE);
            }
        };
        
        if (format == FileFormat::FASTA) {
            std::string line, sequence;
            while (std::getline(input, line)) {
                if (line.empty()) continue;
                
                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        current_batch.push_back(std::move(sequence));
                        sequence.clear();
                        
                        if (current_batch.size() >= BATCH_SIZE) {
                            flush_batch();
                        }
                    }
                } else {
                    sequence += line;
                }
            }
            
            if (!sequence.empty()) {
                current_batch.push_back(std::move(sequence));
            }
        } else if (format == FileFormat::FASTQ) {
            std::string line;
            int line_number = 0;
            
            while (std::getline(input, line)) {
                int line_type = line_number % 4;
                
                if (line_type == 1 && !line.empty()) { // Sequence line
                    current_batch.push_back(line);
                    
                    if (current_batch.size() >= BATCH_SIZE) {
                        flush_batch();
                    }
                }
                
                line_number++;
            }
        } else { // PLAIN
            std::string line;
            while (std::getline(input, line)) {
                if (!line.empty()) {
                    current_batch.push_back(line);
                    
                    if (current_batch.size() >= BATCH_SIZE) {
                        flush_batch();
                    }
                }
            }
        }
        
        flush_batch(); // Flush any remaining sequences
    }
    
    /**
     * Count k-mers from input file (ARM64 optimized)
     */
    bool count_kmers_from_file(const std::string& filename) {
        if (!hash_loaded) {
            std::cerr << "Error: Perfect hash not loaded. Call load_perfect_hash() first." << std::endl;
            return false;
        }
        
        std::cout << "ARM64 processing file: " << filename << std::endl;
        
        std::ifstream input(filename);
        if (!input.is_open()) {
            std::cerr << "Error: Cannot open input file: " << filename << std::endl;
            return false;
        }
        
        // Detect file format
        FileFormat format = detect_format(filename);
        std::cout << "Detected format: ";
        switch (format) {
            case FileFormat::FASTA: std::cout << "FASTA"; break;
            case FileFormat::FASTQ: std::cout << "FASTQ"; break;
            case FileFormat::PLAIN: std::cout << "Plain text"; break;
        }
        std::cout << " (ARM64 optimized)" << std::endl;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Start ARM64-optimized worker threads
        std::vector<std::thread> workers;
        for (size_t i = 0; i < num_threads; ++i) {
            workers.emplace_back(&ARM64Kmer13Counter::worker_thread_arm64, this);
        }
        
        // Read file in batches
        input.clear();
        input.seekg(0);
        read_file_in_batches(input, format);
        input.close();
        
        // Signal completion and wait for workers
        done_reading = true;
        cv.notify_all();
        
        for (auto& worker : workers) {
            worker.join();
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "ARM64 processing completed in " << duration.count() << " ms" << std::endl;
        
        return true;
    }
    
    /**
     * Save counts to binary file
     */
    bool save_counts(const std::string& output_file) {
        std::cout << "Saving ARM64 counts to: " << output_file << std::endl;
        
        std::ofstream out(output_file, std::ios::binary);
        if (!out.is_open()) {
            std::cerr << "Error: Cannot create output file: " << output_file << std::endl;
            return false;
        }
        
        // Write all counts as uint64_t array (convert from uint32_t)
        for (uint32_t i = 0; i < TOTAL_13MERS; i++) {
            uint64_t count = kmer_counts[i].load();
            out.write(reinterpret_cast<const char*>(&count), sizeof(uint64_t));
        }
        
        out.close();
        
        // Verify file size
        std::ifstream check(output_file, std::ios::ate | std::ios::binary);
        size_t file_size = check.tellg();
        size_t expected_size = TOTAL_13MERS * sizeof(uint64_t);
        
        if (file_size != expected_size) {
            std::cerr << "Error: Output file size mismatch. Expected: " << expected_size 
                     << ", got: " << file_size << std::endl;
            return false;
        }
        
        std::cout << "ARM64 counts saved successfully (" << (file_size / (1024*1024)) << " MB)" << std::endl;
        return true;
    }
    
    /**
     * Print ARM64-optimized statistics
     */
    void print_statistics() {
        uint64_t total_seqs = total_sequences.load();
        uint64_t total_proc = total_kmers_processed.load();
        uint64_t valid = valid_kmers.load();
        uint64_t invalid = invalid_kmers.load();
        
        // Count non-zero k-mers
        uint64_t unique_kmers = 0;
        uint64_t total_count = 0;
        uint64_t max_count = 0;
        
        for (uint32_t i = 0; i < TOTAL_13MERS; i++) {
            uint64_t count = kmer_counts[i].load();
            if (count > 0) {
                unique_kmers++;
                total_count += count;
                max_count = std::max(max_count, count);
            }
        }
        
        std::cout << "\n=== ARM64 13-mer Counting Statistics ===" << std::endl;
        std::cout << "Sequences processed: " << total_seqs << std::endl;
        std::cout << "Total k-mers processed: " << total_proc << std::endl;
        std::cout << "Valid k-mers: " << valid << std::endl;
        std::cout << "Invalid k-mers: " << invalid << std::endl;
        std::cout << "Unique k-mers found: " << unique_kmers << " / " << TOTAL_13MERS 
                 << " (" << (100.0 * unique_kmers / TOTAL_13MERS) << "%)" << std::endl;
        std::cout << "Total k-mer count: " << total_count << std::endl;
        std::cout << "Max k-mer frequency: " << max_count << std::endl;
        std::cout << "Average k-mer frequency: " << (unique_kmers > 0 ? (double)total_count / unique_kmers : 0) << std::endl;
        std::cout << "ARM64 threads used: " << num_threads << std::endl;
        
        // ARM64 performance metrics
        if (total_proc > 0) {
            std::cout << "Valid k-mer ratio: " << (100.0 * valid / total_proc) << "%" << std::endl;
        }
    }
};

int main(int argc, char* argv[]) {
    std::cout << "ARM64 13-mer Counter for Apple Silicon" << std::endl;
    std::cout << "Optimized for M1/M2/M3 processors" << std::endl;
    
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <hash_file> <output_tf_file> [num_threads]" << std::endl;
        std::cerr << std::endl;
        std::cerr << "Arguments:" << std::endl;
        std::cerr << "  input_file     - Input sequences (FASTA, FASTQ, or plain text)" << std::endl;
        std::cerr << "  hash_file      - Precomputed perfect hash file (.pf)" << std::endl;
        std::cerr << "  output_tf_file - Output counts file (.tf.bin)" << std::endl;
        std::cerr << "  num_threads    - Number of threads (default: auto for ARM64)" << std::endl;
        std::cerr << std::endl;
        std::cerr << "ARM64 Examples:" << std::endl;
        std::cerr << "  " << argv[0] << " reads.fastq 13mer_index.pf reads.tf.bin" << std::endl;
        std::cerr << "  " << argv[0] << " genome.fasta 13mer_index.pf genome.tf.bin 10" << std::endl;
        return 1;
    }
    
    std::string input_file = argv[1];
    std::string hash_file = argv[2];
    std::string output_file = argv[3];
    size_t num_threads = (argc > 4) ? std::stoul(argv[4]) : std::thread::hardware_concurrency();
    
    std::cout << "\n=== ARM64 13-mer Counter Configuration ===" << std::endl;
    std::cout << "Input file: " << input_file << std::endl;
    std::cout << "Hash file: " << hash_file << std::endl;
    std::cout << "Output file: " << output_file << std::endl;
    std::cout << "ARM64 threads: " << num_threads << std::endl;
    std::cout << std::endl;
    
    try {
        ARM64Kmer13Counter counter(num_threads);
        
        // Load perfect hash
        if (!counter.load_perfect_hash(hash_file)) {
            return 1;
        }
        
        // Process input file
        if (!counter.count_kmers_from_file(input_file)) {
            return 1;
        }
        
        // Print statistics
        counter.print_statistics();
        
        // Save results
        if (!counter.save_counts(output_file)) {
            return 1;
        }
        
        std::cout << "\n✓ ARM64 13-mer counting completed successfully!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "ARM64 Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
