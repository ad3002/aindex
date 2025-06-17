#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>
#include <iomanip>
#include <chrono>
#include "kmers.hpp"

/**
 * Generator of all possible 13-mers
 * Creates a complete set of 4^13 = 67,108,864 k-mers
 */
class All13MerGenerator {
private:
    static const uint32_t TOTAL_13MERS = 67108864; // 4^13
    static const char BASES[4];
    
public:
    /**
     * Generates a k-mer by its numeric index (0 to 4^13-1)
     */
    std::string index_to_kmer(uint32_t index) {
        std::string kmer(13, 'A');
        uint32_t temp = index;
        
        for (int i = 12; i >= 0; i--) {
            kmer[i] = BASES[temp & 3];
            temp >>= 2;
        }
        
        return kmer;
    }
    
    /**
     * Converts a k-mer back to its index for verification
     */
    uint32_t kmer_to_index(const std::string& kmer) {
        return get_dna13_bitset(kmer);
    }
    
    /**
     * Generates all 13-mers and saves them to a file
     * @param output_file - output file name
     * @param with_indices - whether to include numeric indices
     * @param binary_format - whether to save in binary format
     */
    void generate_all_kmers(const std::string& output_file, 
                           bool with_indices = false, 
                           bool binary_format = false) {
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        std::cout << "Generating all " << TOTAL_13MERS << " possible 13-mers..." << std::endl;
        
        if (binary_format) {
            generate_binary(output_file);
        } else {
            generate_text(output_file, with_indices);
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Generation completed in " << duration.count() << " ms" << std::endl;
        std::cout << "Output saved to: " << output_file << std::endl;
    }
    
    /**
     * Generates a text file with k-mers
     */
    void generate_text(const std::string& output_file, bool with_indices) {
        std::ofstream out(output_file);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot create output file: " + output_file);
        }
        
        const uint32_t progress_step = TOTAL_13MERS / 100; // for progress bar
        
        for (uint32_t i = 0; i < TOTAL_13MERS; i++) {
            std::string kmer = index_to_kmer(i);
            
            if (with_indices) {
                out << kmer << "\t" << i << "\n";
            } else {
                out << kmer << "\n";
            }
            
            // Show progress every 1%
            if (i % progress_step == 0) {
                int percent = (i * 100) / TOTAL_13MERS;
                std::cout << "\rProgress: " << percent << "%" << std::flush;
            }
        }
        std::cout << "\rProgress: 100%" << std::endl;
        
        out.close();
    }
    
    /**
     * Generates a binary file with k-mers (more compact)
     */
    void generate_binary(const std::string& output_file) {
        std::ofstream out(output_file + ".bin", std::ios::binary);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot create binary output file");
        }
        
        // Header: magic number, version, k, number of k-mers
        uint32_t magic = 0x4B4D5233; // "KMR3" for 13-mers
        uint32_t version = 1;
        uint32_t k = 13;
        uint32_t count = TOTAL_13MERS;
        
        out.write(reinterpret_cast<const char*>(&magic), sizeof(uint32_t));
        out.write(reinterpret_cast<const char*>(&version), sizeof(uint32_t));
        out.write(reinterpret_cast<const char*>(&k), sizeof(uint32_t));
        out.write(reinterpret_cast<const char*>(&count), sizeof(uint32_t));
        
        const uint32_t progress_step = TOTAL_13MERS / 100;
        
        // Write all k-mers as uint32_t
        for (uint32_t i = 0; i < TOTAL_13MERS; i++) {
            out.write(reinterpret_cast<const char*>(&i), sizeof(uint32_t));
            
            if (i % progress_step == 0) {
                int percent = (i * 100) / TOTAL_13MERS;
                std::cout << "\rProgress: " << percent << "%" << std::flush;
            }
        }
        std::cout << "\rProgress: 100%" << std::endl;
        
        out.close();
    }
    
    /**
     * Checks the correctness of generation (test function)
     */
    void validate_generation(uint32_t num_samples = 1000) {
        std::cout << "Validating generation with " << num_samples << " random samples..." << std::endl;
        
        srand(time(nullptr));
        
        for (uint32_t i = 0; i < num_samples; i++) {
            uint32_t random_index = rand() % TOTAL_13MERS;
            std::string kmer = index_to_kmer(random_index);
            uint32_t back_index = kmer_to_index(kmer);
            
            if (random_index != back_index) {
                std::cerr << "Validation failed at index " << random_index 
                         << ": kmer=" << kmer 
                         << ", back_index=" << back_index << std::endl;
                return;
            }
        }
        
        std::cout << "Validation passed!" << std::endl;
    }
    
    /**
     * Shows statistics about 13-mers
     */
    void print_statistics() {
        std::cout << "\n=== 13-mer Statistics ===" << std::endl;
        std::cout << "Total possible 13-mers: " << TOTAL_13MERS << std::endl;
        std::cout << "Storage as uint32_t: " << (TOTAL_13MERS * sizeof(uint32_t)) / (1024*1024) << " MB" << std::endl;
        std::cout << "Each 13-mer uses: " << sizeof(uint32_t) << " bytes (vs " << 13 << " bytes as string)" << std::endl;
        std::cout << "Compression ratio: " << (13.0 / sizeof(uint32_t)) << "x" << std::endl;
        
        // Show examples
        std::cout << "\nExamples:" << std::endl;
        for (uint32_t i = 0; i < 10; i++) {
            std::string kmer = index_to_kmer(i);
            std::cout << "  Index " << std::setw(8) << i << " -> " << kmer 
                     << " -> " << std::setw(8) << kmer_to_index(kmer) << std::endl;
        }
        
        // Last few
        std::cout << "  ..." << std::endl;
        for (uint32_t i = TOTAL_13MERS - 3; i < TOTAL_13MERS; i++) {
            std::string kmer = index_to_kmer(i);
            std::cout << "  Index " << std::setw(8) << i << " -> " << kmer 
                     << " -> " << std::setw(8) << kmer_to_index(kmer) << std::endl;
        }
    }
};

const char All13MerGenerator::BASES[4] = {'A', 'C', 'G', 'T'};

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            std::cout << "Usage: " << argv[0] << " <output_file> [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  -i, --with-indices    Include numerical indices in output" << std::endl;
            std::cout << "  -b, --binary         Generate binary format" << std::endl;
            std::cout << "  -s, --stats          Show statistics only" << std::endl;
            std::cout << "  -v, --validate       Run validation test" << std::endl;
            std::cout << "\nExamples:" << std::endl;
            std::cout << "  " << argv[0] << " all_13mers.txt" << std::endl;
            std::cout << "  " << argv[0] << " all_13mers.txt -i -b" << std::endl;
            std::cout << "  " << argv[0] << " dummy -s -v" << std::endl;
            return 1;
        }
        
        std::string output_file = argv[1];
        bool with_indices = false;
        bool binary_format = false;
        bool stats_only = false;
        bool validate = false;
        
        // Parse options
        for (int i = 2; i < argc; i++) {
            std::string arg = argv[i];
            if (arg == "-i" || arg == "--with-indices") {
                with_indices = true;
            } else if (arg == "-b" || arg == "--binary") {
                binary_format = true;
            } else if (arg == "-s" || arg == "--stats") {
                stats_only = true;
            } else if (arg == "-v" || arg == "--validate") {
                validate = true;
            }
        }
        
        All13MerGenerator generator;
        
        // Show statistics
        generator.print_statistics();
        
        // Run validation if requested
        if (validate) {
            generator.validate_generation();
        }
        
        // Generate file if not stats only
        if (!stats_only) {
            generator.generate_all_kmers(output_file, with_indices, binary_format);
            
            // Check the size of the created file
            std::ifstream file(output_file, std::ios::ate | std::ios::binary);
            if (file.is_open()) {
                auto size = file.tellg();
                std::cout << "File size: " << (size / (1024*1024)) << " MB" << std::endl;
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
