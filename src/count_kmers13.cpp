//
// Specialized 13-mer counter using precomputed perfect hash
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

#include "emphf/common.hpp"
#include "hash.hpp"
#include "kmers.hpp"

class Kmer13Counter {
private:
    static const uint32_t TOTAL_13MERS = 67108864; // 4^13
    static const int KMER_SIZE = 13;
    
    // Perfect hash components
    HASHER hasher;
    std::atomic<uint64_t>* kmer_counts;
    bool hash_loaded;
    
    // Threading
    std::mutex queue_mutex;
    std::condition_variable cv;
    std::queue<std::string> sequence_queue;
    std::atomic<bool> done_reading{false};
    size_t num_threads;
    
    // Statistics
    std::atomic<uint64_t> total_sequences{0};
    std::atomic<uint64_t> total_kmers_processed{0};
    std::atomic<uint64_t> valid_kmers{0};
    std::atomic<uint64_t> invalid_kmers{0};
    
    emphf::stl_string_adaptor str_adapter;
    
public:
    Kmer13Counter(size_t threads = std::thread::hardware_concurrency()) 
        : kmer_counts(nullptr), hash_loaded(false), num_threads(threads) {
        if (num_threads == 0) num_threads = 1;
        
        // Initialize atomic array for all possible 13-mers
        kmer_counts = new std::atomic<uint64_t>[TOTAL_13MERS];
        for (uint32_t i = 0; i < TOTAL_13MERS; i++) {
            kmer_counts[i] = 0;
        }
    }
    
    ~Kmer13Counter() {
        delete[] kmer_counts;
    }
    
    /**
     * Load precomputed perfect hash
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
            std::cout << "Perfect hash loaded successfully" << std::endl;
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error loading hash: " << e.what() << std::endl;
            return false;
        }
    }
    
    /**
     * Check if k-mer contains only valid nucleotides (ATGC)
     */
    bool is_valid_kmer(const std::string& kmer) {
        for (char c : kmer) {
            if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
                return false;
            }
        }
        return true;
    }
    
    /**
     * Convert string to uppercase and validate
     */
    std::string normalize_sequence(const std::string& seq) {
        std::string normalized;
        normalized.reserve(seq.length());
        
        for (char c : seq) {
            char upper_c = std::toupper(c);
            if (upper_c == 'A' || upper_c == 'T' || upper_c == 'G' || upper_c == 'C') {
                normalized += upper_c;
            } else {
                normalized += 'N'; // Mark invalid nucleotides
            }
        }
        return normalized;
    }
    
    /**
     * Process a single sequence and count its k-mers
     */
    void process_sequence(const std::string& seq) {
        if (seq.length() < KMER_SIZE) return;
        
        std::string normalized = normalize_sequence(seq);
        total_sequences++;
        
        // Count k-mers in this sequence
        for (size_t i = 0; i <= normalized.length() - KMER_SIZE; ++i) {
            std::string kmer = normalized.substr(i, KMER_SIZE);
            total_kmers_processed++;
            
            if (is_valid_kmer(kmer)) {
                // Use perfect hash to get index
                uint64_t hash_idx = hasher.lookup(kmer, str_adapter);
                
                // For 13-mers, we know the hash should map to [0, TOTAL_13MERS)
                if (hash_idx < TOTAL_13MERS) {
                    kmer_counts[hash_idx].fetch_add(1);
                    valid_kmers++;
                } else {
                    std::cerr << "Warning: hash index out of range: " << hash_idx << " for k-mer: " << kmer << std::endl;
                    invalid_kmers++;
                }
            } else {
                invalid_kmers++;
            }
        }
    }
    
    /**
     * Worker thread function
     */
    void worker_thread() {
        while (true) {
            std::unique_lock<std::mutex> lock(queue_mutex);
            cv.wait(lock, [this] { return !sequence_queue.empty() || done_reading; });
            
            if (sequence_queue.empty() && done_reading) {
                break;
            }
            
            if (!sequence_queue.empty()) {
                std::string seq = sequence_queue.front();
                sequence_queue.pop();
                lock.unlock();
                
                process_sequence(seq);
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
     * Read FASTA file
     */
    void read_fasta_file(std::ifstream& input) {
        std::string line, sequence;
        while (std::getline(input, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                if (!sequence.empty()) {
                    std::lock_guard<std::mutex> lock(queue_mutex);
                    sequence_queue.push(sequence);
                    cv.notify_one();
                    sequence.clear();
                }
            } else {
                sequence += line;
            }
        }
        
        // Add last sequence
        if (!sequence.empty()) {
            std::lock_guard<std::mutex> lock(queue_mutex);
            sequence_queue.push(sequence);
            cv.notify_one();
        }
    }
    
    /**
     * Read FASTQ file (SE or PE)
     */
    void read_fastq_file(std::ifstream& input) {
        std::string line;
        int line_number = 0;
        
        while (std::getline(input, line)) {
            int line_type = line_number % 4;
            
            if (line_type == 1) { // Sequence line
                if (!line.empty()) {
                    std::lock_guard<std::mutex> lock(queue_mutex);
                    sequence_queue.push(line);
                    cv.notify_one();
                }
            }
            
            line_number++;
        }
    }
    
    /**
     * Read plain text file (one sequence per line)
     */
    void read_plain_file(std::ifstream& input) {
        std::string line;
        while (std::getline(input, line)) {
            if (!line.empty()) {
                std::lock_guard<std::mutex> lock(queue_mutex);
                sequence_queue.push(line);
                cv.notify_one();
            }
        }
    }
    
    /**
     * Count k-mers from input file
     */
    bool count_kmers_from_file(const std::string& filename) {
        if (!hash_loaded) {
            std::cerr << "Error: Perfect hash not loaded. Call load_perfect_hash() first." << std::endl;
            return false;
        }
        
        std::cout << "Processing file: " << filename << std::endl;
        
        std::ifstream input(filename);
        if (!input.is_open()) {
            std::cerr << "Error: Cannot open input file: " << filename << std::endl;
            return false;
        }
        
        // Detect file format
        FileFormat format = detect_format(filename);
        std::cout << "Detected format: ";
        switch (format) {
            case FileFormat::FASTA: std::cout << "FASTA" << std::endl; break;
            case FileFormat::FASTQ: std::cout << "FASTQ" << std::endl; break;
            case FileFormat::PLAIN: std::cout << "Plain text" << std::endl; break;
        }
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Start worker threads
        std::vector<std::thread> workers;
        for (size_t i = 0; i < num_threads; ++i) {
            workers.emplace_back(&Kmer13Counter::worker_thread, this);
        }
        
        // Read file and feed sequences to workers
        input.clear();
        input.seekg(0);
        
        switch (format) {
            case FileFormat::FASTA:
                read_fasta_file(input);
                break;
            case FileFormat::FASTQ:
                read_fastq_file(input);
                break;
            case FileFormat::PLAIN:
                read_plain_file(input);
                break;
        }
        
        input.close();
        
        // Signal completion and wait for workers
        done_reading = true;
        cv.notify_all();
        
        for (auto& worker : workers) {
            worker.join();
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Processing completed in " << duration.count() << " ms" << std::endl;
        
        return true;
    }
    
    /**
     * Save counts to binary file
     */
    bool save_counts(const std::string& output_file) {
        std::cout << "Saving counts to: " << output_file << std::endl;
        
        std::ofstream out(output_file, std::ios::binary);
        if (!out.is_open()) {
            std::cerr << "Error: Cannot create output file: " << output_file << std::endl;
            return false;
        }
        
        // Write all counts as uint64_t array
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
        
        std::cout << "Counts saved successfully (" << (file_size / (1024*1024)) << " MB)" << std::endl;
        return true;
    }
    
    /**
     * Print statistics
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
        
        std::cout << "\n=== K-mer Counting Statistics ===" << std::endl;
        std::cout << "Sequences processed: " << total_seqs << std::endl;
        std::cout << "Total k-mers processed: " << total_proc << std::endl;
        std::cout << "Valid k-mers: " << valid << std::endl;
        std::cout << "Invalid k-mers: " << invalid << std::endl;
        std::cout << "Unique k-mers found: " << unique_kmers << " / " << TOTAL_13MERS 
                 << " (" << (100.0 * unique_kmers / TOTAL_13MERS) << "%)" << std::endl;
        std::cout << "Total k-mer count: " << total_count << std::endl;
        std::cout << "Max k-mer frequency: " << max_count << std::endl;
        std::cout << "Average k-mer frequency: " << (unique_kmers > 0 ? (double)total_count / unique_kmers : 0) << std::endl;
        std::cout << "Threads used: " << num_threads << std::endl;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <hash_file> <output_tf_file> [num_threads]" << std::endl;
        std::cerr << std::endl;
        std::cerr << "Arguments:" << std::endl;
        std::cerr << "  input_file     - Input sequences (FASTA, FASTQ, or plain text)" << std::endl;
        std::cerr << "  hash_file      - Precomputed perfect hash file (.hash)" << std::endl;
        std::cerr << "  output_tf_file - Output counts file (.tf.bin)" << std::endl;
        std::cerr << "  num_threads    - Number of threads (default: auto)" << std::endl;
        std::cerr << std::endl;
        std::cerr << "Examples:" << std::endl;
        std::cerr << "  " << argv[0] << " reads.fastq 13mer_index.hash reads.tf.bin" << std::endl;
        std::cerr << "  " << argv[0] << " genome.fasta 13mer_index.hash genome.tf.bin 8" << std::endl;
        std::cerr << "  " << argv[0] << " sequences.txt 13mer_index.hash sequences.tf.bin 16" << std::endl;
        return 1;
    }
    
    std::string input_file = argv[1];
    std::string hash_file = argv[2];
    std::string output_file = argv[3];
    size_t num_threads = (argc > 4) ? std::stoul(argv[4]) : std::thread::hardware_concurrency();
    
    std::cout << "=== 13-mer Counter with Perfect Hash ===" << std::endl;
    std::cout << "Input file: " << input_file << std::endl;
    std::cout << "Hash file: " << hash_file << std::endl;
    std::cout << "Output file: " << output_file << std::endl;
    std::cout << "Threads: " << num_threads << std::endl;
    std::cout << std::endl;
    
    try {
        Kmer13Counter counter(num_threads);
        
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
        
        std::cout << "\n✓ 13-mer counting completed successfully!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
