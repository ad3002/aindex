#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/mman.h>
#include <atomic>
#include <mutex>
#include <thread>
#include "emphf/common.hpp"
#include "hash.hpp"
#include <string_view>
#include "helpers.hpp"
#include <fcntl.h>
#include <unistd.h>
#include <unordered_map>
#include <cstring>
#include <string_view>
#include <set>

// Terminology that is used in this file:
//     kmer - std::string
//     ukmer - uint64_t
//     ckmer - char*
//     kid - kmer id, index of kmer in perfect hash
//     pfid - perfect hash id, index of kmer in perfect hash
//     read - sequence of nucleotides from reads file
//     rid - read id is read index in reads file
//     tf - term frequency, number of times kmer appears in reads
//     pos - position in reads file
//     start - start position of read in reads file
//     end - end position of read in reads file
//     local_start - start position of kmer in read

typedef std::atomic<uint8_t> ATOMIC_BOOL;
emphf::stl_string_adaptor str_adapter;

// Define a structure for an interval
struct Interval {
    uint64_t rid;
    uint64_t start;
    uint64_t end;

    bool operator<(const Interval& other) const {
        return start < other.start;
    }
};

// Class to manage intervals
class IntervalTree {
public:
    std::vector<Interval> intervals;

    void addInterval(uint64_t rid, uint64_t start, uint64_t end) {
        intervals.push_back({rid, start, end});
    }

    void sort() {
        std::sort(intervals.begin(), intervals.end());
    }

    std::vector<Interval> query(uint64_t start, uint64_t end) {
        std::vector<Interval> result;
        for (const auto& interval : intervals) {
            if (interval.start <= end && interval.end >= start) {
                result.push_back(interval);
            }
        }
        return result;
    }
};

class UsedReads {
private:
    std::set<uint64_t> read_ids;
    uint64_t read_count = 0;
    uint64_t max_reads = 0;

public:
    UsedReads(uint64_t max_r) : max_reads(max_r) {}

    bool add_read(uint64_t rid) {
        if (read_count >= max_reads) {
            return false;
        }
        
        if (read_ids.find(rid) != read_ids.end()) {
            return true; // Already exists, but we can continue
        }
        
        read_ids.insert(rid);
        read_count++;
        return true;
    }

    bool is_full() const {
        return read_count >= max_reads;
    }

    uint64_t size() const {
        return read_count;
    }

    void clear() {
        read_ids.clear();
        read_count = 0;
    }

    bool contains(uint64_t rid) const {
        return read_ids.find(rid) != read_ids.end();
    }

    std::set<uint64_t> get_reads() const {
        return read_ids;
    }
};

struct Hit {
    uint64_t rid;
    uint64_t start;
    std::string read;
    uint64_t local_pos;
    int ori;
    bool rev;
};

class AindexWrapper {
    uint64_t *positions = nullptr;
    uint64_t *indices = nullptr;
    uint64_t n = 0;
    uint32_t max_tf = 0;
    uint64_t indices_length = 0;

    // 13-mer mode support
    bool is_13mer_mode = false;
    HASHER hasher_13mer;
    uint64_t* tf_array_13mer = nullptr;  // Use uint64_t to match file format
    static const uint32_t TOTAL_13MERS = 67108864; // 4^13
    
    // 13-mer AIndex support (positions)
    uint64_t *positions_13mer = nullptr;
    uint64_t *indices_13mer = nullptr;
    uint64_t n_13mer = 0;
    uint64_t indices_length_13mer = 0;
    
public:
    bool aindex_loaded = false;
    PHASH_MAP *hash_map;
    uint64_t n_reads = 0;
    uint64_t n_kmers = 0;
    
    uint64_t reads_size = 0;
    char *reads = nullptr;
    bool reads_is_mmaped = false;  // Track allocation method

    std::unordered_map<uint64_t, uint32_t> start2rid;
    std::unordered_map<uint64_t, uint64_t> start2end;
    std::vector<uint64_t> start_positions;

    IntervalTree pos_intervalTree;
    
    AindexWrapper() : 
        positions(nullptr),
        indices(nullptr),
        n(0),
        max_tf(0),
        indices_length(0),
        is_13mer_mode(false),
        tf_array_13mer(nullptr),
        positions_13mer(nullptr),
        indices_13mer(nullptr),
        n_13mer(0),
        indices_length_13mer(0),
        aindex_loaded(false),
        hash_map(nullptr),
        n_reads(0),
        n_kmers(0),
        reads_size(0),
        reads(nullptr),
        reads_is_mmaped(false) {}

    ~AindexWrapper() {
        // Safely unmap memory-mapped files
        if (positions != nullptr) {
            munmap(positions, n * sizeof(uint64_t));
            positions = nullptr;
        }
        
        if (indices != nullptr) {
            munmap(indices, indices_length);
            indices = nullptr;
        }
        
        if (reads != nullptr) {
            if (reads_is_mmaped) {
                munmap(reads, reads_size);
            } else {
                delete[] reads;
            }
            reads = nullptr;
        }
        
        if (tf_array_13mer != nullptr) {
            munmap(tf_array_13mer, TOTAL_13MERS * sizeof(uint64_t));
            tf_array_13mer = nullptr;
        }
        
        if (positions_13mer != nullptr) {
            munmap(positions_13mer, n_13mer * sizeof(uint64_t));
            positions_13mer = nullptr;
        }
        
        if (indices_13mer != nullptr) {
            munmap(indices_13mer, indices_length_13mer);
            indices_13mer = nullptr;
        }

        // Safely delete hash_map
        if (hash_map != nullptr) {
            delete hash_map;
            hash_map = nullptr;
        }
    }

    void load(std::string hash_filename, std::string tf_file, std::string kmers_bin_filename, std::string kmers_text_filename){
        // Clean up existing hash_map if it exists
        if (hash_map != nullptr) {
            delete hash_map;
            hash_map = nullptr;
        }
        
        hash_map = new PHASH_MAP();
        // Load perfect hash into hash_map into memory
        emphf::logger() << "Reading index and hash..." << std::endl;
        emphf::logger() << "...files: " << hash_filename << std::endl;
        emphf::logger() << "...files: " << tf_file << std::endl;
        emphf::logger() << "...files: " << kmers_bin_filename << std::endl;
        emphf::logger() << "...files: " << kmers_text_filename << std::endl;
        load_hash(*hash_map, hash_filename, tf_file, kmers_bin_filename, kmers_text_filename);
        n_kmers = hash_map->n;
        emphf::logger() << "\tDone" << std::endl;
    }

    void load_hash_file(std::string hash_filename, std::string tf_file, std::string kmers_bin_filename, std::string kmers_text_filename) {
        emphf::logger() << "Loading hash with all files..." << std::endl;
        
        // Clean up existing hash_map if it exists
        if (hash_map != nullptr) {
            delete hash_map;
            hash_map = nullptr;
        }
        
        hash_map = new PHASH_MAP();
        load_hash(*hash_map, hash_filename, tf_file, kmers_bin_filename, kmers_text_filename);
        n_kmers = hash_map->n;
    }

    void load_reads_index(const std::string& index_file) {
        std::ifstream fin(index_file, std::ios::in);
        if (!fin.is_open()) {
            std::cerr << "Error opening index file: " << index_file << std::endl;
            std::terminate();
        }

        n_reads = 0;
        uint64_t rid, start_pos, end_pos;
        while (fin >> rid >> start_pos >> end_pos) {
            pos_intervalTree.addInterval(rid, start_pos, end_pos+1);
            start2rid[start_pos] = rid;
            start_positions.push_back(start_pos);
            start2end[start_pos] = end_pos;
            n_reads++;
        }

        fin.close();
    }

    void load_reads(std::string reads_file) {
        // Clean up existing reads mapping if it exists
        if (reads != nullptr) {
            if (reads_is_mmaped) {
                munmap(reads, reads_size);
            } else {
                delete[] reads;
            }
            reads = nullptr;
            reads_size = 0;
        }
        
        // Memory map reads
        emphf::logger() << "Memory mapping reads file..." << std::endl;
        std::ifstream fout(reads_file, std::ios::in | std::ios::binary);
        fout.seekg(0, std::ios::end);
        uint64_t length = fout.tellg();
        fout.close();

        FILE* in = std::fopen(reads_file.c_str(), "rb");
        if (in == nullptr) {
            std::cerr << "Failed to open reads file: " << reads_file << std::endl;
            return;
        }
        
        reads = (char*)mmap(NULL, length, PROT_READ|PROT_WRITE, MAP_PRIVATE, fileno(in), 0);
        if (reads == MAP_FAILED) {
            std::cerr << "Failed to mmap reads file" << std::endl;
            reads = nullptr;
            fclose(in);
            return;
        }
        fclose(in);

        reads_size = length;
        reads_is_mmaped = true;

        emphf::logger() << "\tbuilding start pos index over reads: " << std::endl;
        std::string index_file = reads_file.substr(0, reads_file.find_last_of(".")) + ".ridx";
        load_reads_index(index_file);
        emphf::logger() << "\tDone" << std::endl;
    }

    void load_reads_in_memory(std::string reads_file) {
        // Clean up existing reads if it exists
        if (reads != nullptr) {
            if (reads_is_mmaped) {
                munmap(reads, reads_size);
            } else {
                delete[] reads;
            }
            reads = nullptr;
            reads_size = 0;
        }
        
        // Load reads into memory
        emphf::logger() << "Loading reads file into memory..." << std::endl;
        std::ifstream fin(reads_file, std::ios::in | std::ios::binary);
        if (!fin) {
            std::cerr << "Failed to open file" << std::endl;
            return;
        }

        fin.seekg(0, std::ios::end);
        uint64_t length = fin.tellg();
        fin.seekg(0, std::ios::beg);

        reads = new char[length];
        fin.read(reads, length);
        fin.close();

        reads_size = length;
        reads_is_mmaped = false;  // This is allocated with new[]

        emphf::logger() << "\tbuilding start pos index over reads: " << std::endl;
        std::string index_file = reads_file.substr(0, reads_file.find_last_of(".")) + ".ridx";
        load_reads_index(index_file);
        emphf::logger() << "\tDone" << std::endl;
    }

    void load_aindex(std::string index_file, std::string indices_file, uint32_t _max_tf) {
        // Load aindex.

        n = hash_map->n;
        max_tf = _max_tf;

        emphf::logger() << "Reading aindex.indices.bin array..." << std::endl;

        std::ifstream fin_temp(indices_file, std::ios::in | std::ios::binary);
        fin_temp.seekg(0, std::ios::end);
        uint64_t length = fin_temp.tellg();
        fin_temp.close();

        FILE* in1 = std::fopen(indices_file.c_str(), "rb");
        indices = (uint64_t*)mmap(NULL, length, PROT_READ|PROT_WRITE, MAP_PRIVATE, fileno(in1), 0);
        if (indices == nullptr) {
            std::cerr << "Failed position loading" << std::endl;
            exit(10);
        }
        fclose(in1);
        indices_length = length;
        emphf::logger() << "\tindices length: " << indices_length << std::endl;
        emphf::logger() << "\tDone" << std::endl;

        emphf::logger() << "Reading aindex.index.bin array..." << std::endl;

        std::ifstream fout6(index_file, std::ios::in | std::ios::binary);
        fout6.seekg(0, std::ios::end);
        length = fout6.tellg();
        fout6.close();

        emphf::logger() << "\tpositions length: " << length << std::endl;
        FILE* in = std::fopen(index_file.c_str(), "rb");
        positions = (uint64_t*)mmap(NULL, length, PROT_READ|PROT_WRITE, MAP_PRIVATE, fileno(in), 0);
        if (positions == nullptr) {
            std::cerr << "Failed position loading" << std::endl;
            exit(10);
        }
        fclose(in);
        this->aindex_loaded = true;
        emphf::logger() << "\tDone" << std::endl;
    }

    void load_13mer_index(const std::string& hash_file, const std::string& tf_file) {
        emphf::logger() << "Loading 13-mer index..." << std::endl;
        emphf::logger() << "Hash file: " << hash_file << std::endl;
        emphf::logger() << "TF file: " << tf_file << std::endl;
        
        // Load hash using ifstream
        std::ifstream in(hash_file, std::ios::binary);
        if (!in) {
            std::cerr << "Failed to open hash file: " << hash_file << std::endl;
            std::terminate();
        }
        hasher_13mer.load(in);
        in.close();
        
        // Memory map tf array
        FILE* tf_in = std::fopen(tf_file.c_str(), "rb");
        if (!tf_in) {
            std::cerr << "Failed to open tf file: " << tf_file << std::endl;
            std::terminate();
        }
        
        tf_array_13mer = (uint64_t*)mmap(NULL, TOTAL_13MERS * sizeof(uint64_t), 
                                        PROT_READ, MAP_SHARED, fileno(tf_in), 0);
        if (tf_array_13mer == MAP_FAILED) {
            std::cerr << "Failed to mmap tf file" << std::endl;
            std::terminate();
        }
        fclose(tf_in);
        
        is_13mer_mode = true;
        n_kmers = TOTAL_13MERS;
        
        emphf::logger() << "13-mer index loaded successfully" << std::endl;
    }
    
    void load_13mer_aindex(const std::string& index_file, const std::string& indices_file) {
        emphf::logger() << "Loading 13-mer AIndex files..." << std::endl;
        emphf::logger() << "Index file: " << index_file << std::endl;
        emphf::logger() << "Indices file: " << indices_file << std::endl;
        
        
        // Load indices.bin file
        std::ifstream indices_in(indices_file, std::ios::in | std::ios::binary);
        if (!indices_in) {
            std::cerr << "Failed to open indices file: " << indices_file << std::endl;
            std::terminate();
        }
        indices_in.seekg(0, std::ios::end);
        indices_length_13mer = indices_in.tellg();
        indices_in.close();
        
        emphf::logger() << "\tIndices length: " << indices_length_13mer << std::endl;
        
        FILE* indices_fp = std::fopen(indices_file.c_str(), "rb");
        if (!indices_fp) {
            std::cerr << "Failed to open indices file for mmap: " << indices_file << std::endl;
            std::terminate();
        }
        
        indices_13mer = (uint64_t*)mmap(NULL, indices_length_13mer, PROT_READ, MAP_SHARED, fileno(indices_fp), 0);
        if (indices_13mer == MAP_FAILED) {
            std::cerr << "Failed to mmap indices file" << std::endl;
            std::terminate();
        }
        fclose(indices_fp);
        this->aindex_loaded = true;
        emphf::logger() << "13-mer AIndex loaded successfully" << std::endl;
    }
    
    // String overloads that delegate to string_view versions
    bool is_13mer(const std::string& kmer) const {
        return is_13mer(std::string_view(kmer));
    }
    
    bool is_23mer(const std::string& kmer) const {
        return is_23mer(std::string_view(kmer));
    }
    
    uint32_t get_tf_value_13mer(const std::string& kmer) {
        if (kmer.length() != 13) {
            return 0;
        }
        
        // Validate k-mer (only ATGC)
        for (char c : kmer) {
            if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
                return 0;
            }
        }
        
        emphf::stl_string_adaptor str_adapter;
        
        // Use perfect hash lookup
        uint64_t hash_id = hasher_13mer.lookup(kmer, str_adapter);
        if (hash_id < TOTAL_13MERS) {
            return tf_array_13mer[hash_id];
        }
        
        return 0;
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

    /**
     * Get total TF value for 13-mer (forward + reverse complement)
     */
    uint64_t get_total_tf_value_13mer(const std::string& kmer) {
        if (!is_13mer_mode) {
            std::cerr << "Error: 13-mer mode not enabled" << std::endl;
            return 0;
        }
        
        if (kmer.length() != 13) {
            std::cerr << "Error: k-mer length must be 13, got: " << kmer.length() << std::endl;
            return 0;
        }
        
        // Get TF for forward k-mer
        uint64_t hash_idx = hasher_13mer.lookup(kmer, str_adapter);
        uint64_t tf_forward = tf_array_13mer[hash_idx];
        
        // Get TF for reverse complement
        std::string rc_kmer = get_reverse_complement_13mer(kmer);
        uint64_t rc_hash_idx = hasher_13mer.lookup(rc_kmer, str_adapter);
        uint64_t tf_reverse = tf_array_13mer[rc_hash_idx];
        
        return tf_forward + tf_reverse;
    }

    /**
     * Get total TF values for multiple 13-mers (forward + reverse complement)
     */
    std::vector<uint64_t> get_total_tf_values_13mer(const std::vector<std::string>& kmers) {
        if (!is_13mer_mode) {
            std::cerr << "Error: 13-mer mode not enabled" << std::endl;
            return std::vector<uint64_t>(kmers.size(), 0);
        }
        
        std::vector<uint64_t> total_tfs;
        total_tfs.reserve(kmers.size());
        
        for (const auto& kmer : kmers) {
            total_tfs.push_back(get_total_tf_value_13mer(kmer));
        }
        
        return total_tfs;
    }

    /**
     * Get TF values for 13-mer in both directions (forward, reverse complement)
     */
    std::pair<uint64_t, uint64_t> get_tf_both_directions_13mer(const std::string& kmer) {
        if (!is_13mer_mode) {
            std::cerr << "Error: 13-mer mode not enabled" << std::endl;
            return std::make_pair(0, 0);
        }
        
        if (kmer.length() != 13) {
            std::cerr << "Error: k-mer length must be 13, got: " << kmer.length() << std::endl;
            return std::make_pair(0, 0);
        }
        
        // Get TF for forward k-mer
        uint64_t hash_idx = hasher_13mer.lookup(kmer, str_adapter);
        uint64_t tf_forward = tf_array_13mer[hash_idx];
        
        // Get TF for reverse complement
        std::string rc_kmer = get_reverse_complement_13mer(kmer);
        uint64_t rc_hash_idx = hasher_13mer.lookup(rc_kmer, str_adapter);
        uint64_t tf_reverse = tf_array_13mer[rc_hash_idx];
        
        return std::make_pair(tf_forward, tf_reverse);
    }

    /**
     * Get TF values for multiple 13-mers in both directions
     * Returns vector of pairs: [(forward_tf, reverse_tf), ...]
     */
    std::vector<std::pair<uint64_t, uint64_t>> get_tf_both_directions_13mer_batch(const std::vector<std::string>& kmers) {
        if (!is_13mer_mode) {
            std::cerr << "Error: 13-mer mode not enabled" << std::endl;
            return std::vector<std::pair<uint64_t, uint64_t>>(kmers.size(), std::make_pair(0, 0));
        }
        
        std::vector<std::pair<uint64_t, uint64_t>> results;
        results.reserve(kmers.size());
        
        for (const auto& kmer : kmers) {
            results.push_back(get_tf_both_directions_13mer(kmer));
        }
        
        return results;
    }

    uint32_t get_tf_value_23mer(const std::string& kmer) {
        uint64_t ukmer = get_dna23_bitset(kmer);
        auto h1 = hash_map->hasher.lookup(kmer, str_adapter);
        if (h1 >= hash_map->n || hash_map->checker[h1] != ukmer) {
            std::string rev_kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
            uint64_t urev_kmer = reverseDNA(ukmer);
            get_bitset_dna23(urev_kmer, rev_kmer);
            auto h2 = hash_map->hasher.lookup(rev_kmer, str_adapter);
            if (h2 >= hash_map->n || hash_map->checker[h2] != urev_kmer) {
                return 0;
            } else {
                return hash_map->tf_values[h2];
            }
        } else {
            return hash_map->tf_values[h1];
        }
        return 0;
    }

    std::vector<uint64_t> get_hash_values(std::vector<std::string> kmers) {
        std::vector<uint64_t> hash_values;
        for (const auto& kmer : kmers) {
            uint64_t kmer_id = hash_map->hasher.lookup(kmer, str_adapter);
            hash_values.push_back(kmer_id);
        }
        return hash_values;
    }

    uint64_t get_hash_value(std::string kmer) {
        uint64_t kmer_id = hash_map->hasher.lookup(kmer, str_adapter);
        return kmer_id;
    }

    // General get_tf_value method that auto-detects mode
    uint32_t get_tf_value(const std::string& kmer) {
        if (is_13mer_mode) {
            return get_tf_value_13mer(kmer);
        } else {
            return get_tf_value_23mer(kmer);
        }
    }

    // General get_tf_values method that auto-detects mode
    std::vector<uint32_t> get_tf_values(const std::vector<std::string>& kmers) {
        if (is_13mer_mode) {
            return get_tf_values_13mer(kmers);
        } else {
            std::vector<uint32_t> tf_values;
            tf_values.reserve(kmers.size());
            for (const auto& kmer : kmers) {
                tf_values.push_back(get_tf_value_23mer(kmer));
            }
            return tf_values;
        }
    }

    std::string get_read_by_rid(uint64_t rid) {
        if (start_positions.size() <= rid) {
            return "";
        }
        uint64_t start = start_positions[rid];
        uint64_t end = start2end[start];
        
        std::string read(reads + start, end - start);
        return read;
    }

    std::string get_read(uint64_t start, uint64_t end, bool revcomp = false) {
        if (start >= reads_size || end >= reads_size || start > end) {
            return "";
        }
        
        std::string read(reads + start, end - start);
        if (revcomp) {
            std::string rev_read = "";
            for (int i = read.length() - 1; i >= 0; i--) {
                switch (read[i]) {
                    case 'A': rev_read += 'T'; break;
                    case 'T': rev_read += 'A'; break;
                    case 'C': rev_read += 'G'; break;
                    case 'G': rev_read += 'C'; break;
                    case 'N': rev_read += 'N'; break;
                    default: rev_read += read[i]; break;
                }
            }
            return rev_read;
        }
        return read;
    }

    uint64_t get_kid_by_kmer(std::string kmer) {
        uint64_t ukmer = get_dna23_bitset(kmer);
        auto h1 = hash_map->hasher.lookup(kmer, str_adapter);
        if (h1 >= hash_map->n || hash_map->checker[h1] != ukmer) {
            std::string rev_kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
            uint64_t urev_kmer = reverseDNA(ukmer);
            get_bitset_dna23(urev_kmer, rev_kmer);
            auto h2 = hash_map->hasher.lookup(rev_kmer, str_adapter);
            if (h2 >= hash_map->n || hash_map->checker[h2] != urev_kmer) {
                return 0;
            } else {
                return h2;
            }
        } else {
            return h1;
        }
    }

    std::string get_kmer_by_kid(uint64_t kid) {
        if (kid >= hash_map->n) {
            return "";
        }
        uint64_t ukmer = hash_map->checker[kid];
        return get_bitset_dna23(ukmer);
    }

    uint64_t get_strand(std::string kmer) {
        uint64_t ukmer = get_dna23_bitset(kmer);
        auto h1 = hash_map->hasher.lookup(kmer, str_adapter);
        if (h1 >= hash_map->n || hash_map->checker[h1] != ukmer) {
            std::string rev_kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
            uint64_t urev_kmer = reverseDNA(ukmer);
            get_bitset_dna23(urev_kmer, rev_kmer);
            auto h2 = hash_map->hasher.lookup(rev_kmer, str_adapter);
            if (h2 >= hash_map->n || hash_map->checker[h2] != urev_kmer) {
                return 0; // NOT_FOUND
            } else {
                return 2; // REVERSE
            }
        } else {
            return 1; // FORWARD
        }
    }

    std::tuple<uint64_t, std::string, std::string> get_kmer_info(uint64_t kid) {
        if (kid >= hash_map->n) {
            return std::make_tuple(0, "", "");
        }
        uint64_t ukmer = hash_map->checker[kid];
        std::string kmer = get_bitset_dna23(ukmer);
        uint64_t urev_kmer = reverseDNA(ukmer);
        std::string rkmer = get_bitset_dna23(urev_kmer);
        // Use load() to get a non-atomic value from the atomic tf_values
        uint64_t tf_value = static_cast<uint64_t>(hash_map->tf_values[kid].load());
        return std::make_tuple(tf_value, kmer, rkmer);
    }

    uint64_t get_rid(uint64_t pos) {
        if (!aindex_loaded || pos_intervalTree.intervals.empty()) {
            return 0;
        }
        
        try {
            std::vector<Interval> overlapping = pos_intervalTree.query(pos, pos + 1);
            if (!overlapping.empty()) {
                return overlapping[0].rid;
            }
        } catch (...) {
            // Handle any exceptions in interval tree query
            return 0;
        }
        return 0;
    }

    uint64_t get_start(uint64_t pos) {
        if (!aindex_loaded || pos_intervalTree.intervals.empty()) {
            return 0;
        }
        
        try {
            std::vector<Interval> overlapping = pos_intervalTree.query(pos, pos + 1);
            if (!overlapping.empty()) {
                return overlapping[0].start;
            }
        } catch (...) {
            // Handle any exceptions in interval tree query
            return 0;
        }
        return 0;
    }

    inline bool is_13mer(std::string_view kmer) const noexcept {
        return kmer.length() == 13;
    }
    
    inline bool is_23mer(std::string_view kmer) const noexcept {
        return kmer.length() == 23;
    }

    // ------------- fast version for 23-mers -----------------------
    inline std::vector<uint64_t>
    get_positions_23mer(std::string_view kmer) const noexcept
    {
        std::vector<uint64_t> out;

        if (!hash_map || !positions || !indices)             // early-exit
            return out;

        const uint64_t h1 = hash_map->get_pfid(kmer);
        if (h1 + 1 >= indices_length)                        // bad bucket
            return out;

        const uint64_t* beg = positions + indices[h1];
        const uint64_t* end = positions + indices[h1 + 1];

        out.reserve(end - beg);                              // avoid realloc

        for (const uint64_t* p = beg; p != end; ++p)
            if (uint64_t pos = *p; pos)                      // skip zeros
                out.emplace_back(pos - 1);

        return out;                                          // NRVO / RVO
    }

    // ------------- public dispatcher ------------------------------
    inline std::vector<uint64_t>
    get_positions(std::string_view kmer) const noexcept
    {
        return is_13mer(kmer) ? get_positions_13mer(kmer)
            : is_23mer(kmer) ? get_positions_23mer(kmer)
                            : std::vector<uint64_t>{};
    }

    // ------------- string overloads for pybind11 compatibility ----
    std::vector<uint64_t> get_positions(const std::string& kmer) {
        return get_positions(std::string_view(kmer));
    }
    
    std::vector<uint64_t> get_positions_13mer(const std::string& kmer) {
        return get_positions_13mer(std::string_view(kmer));
    }
    
    std::vector<uint64_t> get_positions_23mer(const std::string& kmer) {
        return get_positions_23mer(std::string_view(kmer));
    }

    uint64_t get_hash_size() {
        if (is_13mer_mode) {
            return TOTAL_13MERS;
        }
        return hash_map ? hash_map->n : 0;
    }

    uint64_t get_reads_size() {
        return reads_size;
    }

    void check_get_reads_se_by_kmer(uint64_t kmer_id, UsedReads& used_reads, std::vector<Hit>& hits) {
        if (!aindex_loaded) {
            emphf::logger() << "Aindex not loaded!" << std::endl;
            return;
        }

        uint64_t start_pos = positions[kmer_id];
        uint64_t end_pos = (kmer_id == n - 1) ? indices_length : positions[kmer_id + 1];
        
        for (uint64_t i = start_pos; i < end_pos; i++) {
            uint64_t pos = indices[i];
            
            std::vector<Interval> overlapping = pos_intervalTree.query(pos, pos + Settings::K - 1);
            
            for (const auto& interval : overlapping) {
                if (used_reads.is_full()) {
                    return;
                }
                
                if (!used_reads.add_read(interval.rid)) {
                    continue;
                }
                
                std::string read = get_read_by_rid(interval.rid);
                if (read.empty()) continue;
                
                uint64_t local_pos = pos - interval.start;
                if (local_pos + Settings::K <= read.length()) {
                    Hit hit;
                    hit.rid = interval.rid;
                    hit.start = interval.start;
                    hit.read = read;
                    hit.local_pos = local_pos;
                    hit.ori = 1;
                    hit.rev = false;
                    hits.push_back(hit);
                }
            }
        }
    }

    std::vector<std::string> get_reads_se_by_kmer(std::string kmer, uint64_t max_reads) {
        std::vector<std::string> result;
        UsedReads used_reads(max_reads);
        std::vector<Hit> hits;
        
        uint64_t kmer_id = hash_map->hasher.lookup(kmer, str_adapter);
        check_get_reads_se_by_kmer(kmer_id, used_reads, hits);
        
        for (const auto& hit : hits) {
            result.push_back(hit.read);
        }
        
        return result;
    }

    void debug_kmer_tf_values() {
        std::vector<uint64_t> h1_values = {1, 10, 100, 1000, 10000, 100000};
        UsedReads used_reads(100);
        std::vector<Hit> hits;

        for (uint64_t h1: h1_values) {
            if (h1 >= n_kmers) continue;

            uint64_t h1_kmer = hash_map->checker[h1];
            std::string kmer = get_bitset_dna23(h1_kmer);
            hits.clear();
            check_get_reads_se_by_kmer(h1, used_reads, hits);

            uint64_t max_pos = 0;

            for (auto hit: hits) {
                max_pos = std::max(max_pos, hit.local_pos);
                std::string subkmer = hit.read.substr(hit.local_pos, Settings::K);
                assert(subkmer == kmer);
                std::cout << kmer << " " << subkmer << " " << h1 << " " << hash_map->tf_values[h1] << std::endl;
            }
        }
    }
    
    // Additional 13-mer functions
    std::vector<uint32_t> get_tf_values_13mer(const std::vector<std::string>& kmers) {
        std::vector<uint32_t> tf_values;
        tf_values.reserve(kmers.size());
        
        if (!is_13mer_mode) {
            // Fill with zeros if not in 13-mer mode
            tf_values.resize(kmers.size(), 0);
            return tf_values;
        }
        
        emphf::stl_string_adaptor str_adapter;
        
        for (const auto& kmer : kmers) {
            if (!is_13mer(std::string_view(kmer))) {
                tf_values.push_back(0);
                continue;
            }
            
            // Validate k-mer (only ATGC)
            bool valid = true;
            for (char c : kmer) {
                if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
                    valid = false;
                    break;
                }
            }
            
            if (!valid) {
                tf_values.push_back(0);
                continue;
            }
            
            // Use perfect hash lookup
            uint64_t hash_id = hasher_13mer.lookup(kmer, str_adapter);
            if (hash_id < TOTAL_13MERS) {
                tf_values.push_back(tf_array_13mer[hash_id]);
            } else {
                tf_values.push_back(0);
            }
        }
        
        return tf_values;
    }
    
    // Direct access to 13-mer tf array
    std::vector<uint32_t> get_13mer_tf_array() {
        if (!is_13mer_mode) {
            return std::vector<uint32_t>();
        }
        
        std::vector<uint32_t> result(tf_array_13mer, tf_array_13mer + TOTAL_13MERS);
        return result;
    }
    
    // Get tf value by direct array index (for 13-mers)
    uint32_t get_tf_by_index_13mer(uint64_t index) {
        if (!is_13mer_mode || index >= TOTAL_13MERS) {
            return 0;
        }
        return tf_array_13mer[index];
    }
    
    // Get statistics about loaded index
    std::string get_index_info() {
        std::string info = "Index Info:\n";
        if (is_13mer_mode && tf_array_13mer != nullptr) {
            info += "Mode: 13-mer\n";
            info += "Total k-mers: " + std::to_string(TOTAL_13MERS) + "\n";
            
            // Count non-zero entries (with safety check)
            uint64_t non_zero_count = 0;
            uint64_t total_count = 0;
            for (uint64_t i = 0; i < TOTAL_13MERS; i++) {
                if (tf_array_13mer[i] > 0) {
                    non_zero_count++;
                    total_count += tf_array_13mer[i];  // uint64_t
                }
            }
            info += "Non-zero entries: " + std::to_string(non_zero_count) + "\n";
            info += "Total k-mer count: " + std::to_string(total_count) + "\n";
        } else if (hash_map != nullptr) {
            info += "Mode: 23-mer\n";
            info += "Total k-mers: " + std::to_string(hash_map->n) + "\n";
        } else {
            info += "Mode: No index loaded\n";
        }
        
        if (aindex_loaded) {
            info += "AIndex: Loaded\n";
            info += "Reads: " + std::to_string(n_reads) + "\n";
        } else {
            info += "AIndex: Not loaded\n";
        }
        
        return info;
    }
    
    /**
     * Get statistics about the 13-mer index
     */
    std::map<std::string, uint64_t> get_13mer_statistics() {
        std::map<std::string, uint64_t> stats;
        
        if (!is_13mer_mode) {
            std::cerr << "Error: 13-mer mode not enabled" << std::endl;
            return stats;
        }
        
        uint64_t total_kmers = TOTAL_13MERS;
        uint64_t non_zero_kmers = 0;
        uint64_t max_frequency = 0;
        uint64_t total_count = 0;
        
        for (uint32_t i = 0; i < TOTAL_13MERS; i++) {
            uint64_t tf = tf_array_13mer[i];
            if (tf > 0) {
                non_zero_kmers++;
                total_count += tf;
                if (tf > max_frequency) {
                    max_frequency = tf;
                }
            }
        }
        
        stats["total_kmers"] = total_kmers;
        stats["non_zero_kmers"] = non_zero_kmers;
        stats["max_frequency"] = max_frequency;
        stats["total_count"] = total_count;
        
        return stats;
    }
    
    inline std::vector<uint64_t> get_positions_13mer(std::string_view kmer) const {
        std::vector<uint64_t> result;
        
        if (!is_13mer_mode || !is_13mer(kmer) || positions_13mer == nullptr || indices_13mer == nullptr) {
            return result;
        }
        
        // Validate k-mer (only ATGC)
        for (char c : kmer) {
            if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
                return result;
            }
        }
        
        emphf::stl_string_adaptor str_adapter;
        
        // Use perfect hash lookup
        uint64_t hash_id = hasher_13mer.lookup(kmer, str_adapter);
        if (hash_id < TOTAL_13MERS) {
            // Get positions for this k-mer
            uint64_t start_idx = indices_13mer[hash_id];
            uint64_t end_idx = indices_13mer[hash_id + 1];
            
            for (uint64_t i = start_idx; i < end_idx && i < n_13mer; ++i) {
                if (positions_13mer[i] > 0) {
                    result.push_back(positions_13mer[i] - 1); // Convert from 1-based to 0-based
                }
            }
        }
        
        return result;
    }
    
    void load_from_prefix_23mer(const std::string& prefix, const std::string& reads_file = "") {
        emphf::logger() << "Loading 23-mer index from prefix: " << prefix << std::endl;
        
        // Construct file paths
        std::string pf_file = prefix + ".pf";
        std::string tf_file = prefix + ".tf.bin";
        std::string kmers_bin_file = prefix + ".kmers.bin";
        std::string kmers_text_file = prefix + ".txt";
        
        // Check required files exist
        std::vector<std::string> required_files = {pf_file, tf_file, kmers_bin_file};
        for (const auto& file : required_files) {
            std::ifstream test(file);
            if (!test.good()) {
                std::cerr << "Required file not found: " << file << std::endl;
                std::terminate();
            }
        }
        
        // Load hash
        emphf::logger() << "Loading 23-mer hash..." << std::endl;
        load_hash_file(pf_file, tf_file, kmers_bin_file, kmers_text_file);
        emphf::logger() << "23-mer hash loaded successfully" << std::endl;
        
        // Load reads if file provided
        if (!reads_file.empty()) {
            load_reads(reads_file);
            emphf::logger() << "Reads loaded from: " << reads_file << std::endl;
        }
    }
    
    void load_aindex_from_prefix_23mer(const std::string& prefix, uint32_t max_tf, const std::string& reads_file = "") {
        emphf::logger() << "Loading 23-mer AIndex from prefix: " << prefix << std::endl;
        
        // Construct file paths
        std::string index_file = prefix + ".index.bin";
        std::string indices_file = prefix + ".indices.bin";
        
        // Check required files exist
        std::vector<std::string> required_files = {index_file, indices_file};
        for (const auto& file : required_files) {
            std::ifstream test(file);
            if (!test.good()) {
                std::cerr << "Required AIndex file not found: " << file << std::endl;
                std::terminate();
            }
        }
        
        // Load aindex
        load_aindex(index_file, indices_file, max_tf);
        emphf::logger() << "23-mer AIndex loaded successfully" << std::endl;
        
        // Load reads if file provided and not already loaded
        if (!reads_file.empty() && reads == nullptr) {
            load_reads(reads_file);
            emphf::logger() << "Reads loaded from: " << reads_file << std::endl;
        }
    }
    
    void load_from_prefix_13mer(const std::string& prefix, const std::string& reads_file = "") {
        emphf::logger() << "Loading 13-mer index from prefix: " << prefix << std::endl;
        
        // Construct file paths
        std::string hash_file = prefix + ".pf";
        std::string tf_file = prefix + ".tf.bin";
        
        // Check required files exist
        std::vector<std::string> required_files = {hash_file, tf_file};
        for (const auto& file : required_files) {
            std::ifstream test(file);
            if (!test.good()) {
                std::cerr << "Required 13-mer file not found: " << file << std::endl;
                std::terminate();
            }
        }
        
        // Load 13-mer index
        load_13mer_index(hash_file, tf_file);
        emphf::logger() << "13-mer index loaded successfully" << std::endl;
        
        // Load reads if file provided
        if (!reads_file.empty()) {
            load_reads(reads_file);
            emphf::logger() << "Reads loaded from: " << reads_file << std::endl;
        }
    }
    
    void load_aindex_from_prefix_13mer(const std::string& prefix, const std::string& reads_file = "") {
        emphf::logger() << "Loading 13-mer AIndex from prefix: " << prefix << std::endl;
        
        // Construct file paths
        std::string index_file = prefix + ".index.bin";
        std::string indices_file = prefix + ".indices.bin";
        
        // Check required files exist
        std::vector<std::string> required_files = {index_file, indices_file};
        for (const auto& file : required_files) {
            std::ifstream test(file);
            if (!test.good()) {
                std::cerr << "Required 13-mer AIndex file not found: " << file << std::endl;
                std::terminate();
            }
        }
        
        // Load 13-mer aindex
        load_13mer_aindex(index_file, indices_file);
        emphf::logger() << "13-mer AIndex loaded successfully" << std::endl;
        
        // Load reads if file provided and not already loaded
        if (!reads_file.empty() && reads == nullptr) {
            load_reads(reads_file);
            emphf::logger() << "Reads loaded from: " << reads_file << std::endl;
        }
    }
    
    // 23-mer specific methods
    std::vector<uint32_t> get_tf_values_23mer(const std::vector<std::string>& kmers) {
        std::vector<uint32_t> tf_values;
        tf_values.reserve(kmers.size());
        
        for (const auto& kmer : kmers) {
            tf_values.push_back(get_tf_value_23mer(kmer));
        }
        
        return tf_values;
    }
    
    uint64_t get_total_tf_value_23mer(const std::string& kmer) {
        if (kmer.length() != 23) {
            return 0;
        }
        
        uint32_t forward_tf = get_tf_value_23mer(kmer);
        
        // Get reverse complement
        std::string rev_kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
        uint64_t ukmer = get_dna23_bitset(kmer);
        uint64_t urev_kmer = reverseDNA(ukmer);
        get_bitset_dna23(urev_kmer, rev_kmer);
        
        uint32_t reverse_tf = get_tf_value_23mer(rev_kmer);
        
        return static_cast<uint64_t>(forward_tf) + static_cast<uint64_t>(reverse_tf);
    }
    
    std::vector<uint64_t> get_total_tf_values_23mer(const std::vector<std::string>& kmers) {
        std::vector<uint64_t> total_tfs;
        total_tfs.reserve(kmers.size());
        
        for (const auto& kmer : kmers) {
            total_tfs.push_back(get_total_tf_value_23mer(kmer));
        }
        
        return total_tfs;
    }
    
    std::pair<uint32_t, uint32_t> get_tf_both_directions_23mer(const std::string& kmer) {
        if (kmer.length() != 23) {
            return {0, 0};
        }
        
        uint32_t forward_tf = get_tf_value_23mer(kmer);
        
        // Get reverse complement
        std::string rev_kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
        uint64_t ukmer = get_dna23_bitset(kmer);
        uint64_t urev_kmer = reverseDNA(ukmer);
        get_bitset_dna23(urev_kmer, rev_kmer);
        
        uint32_t reverse_tf = get_tf_value_23mer(rev_kmer);
        
        return {forward_tf, reverse_tf};
    }
    
    std::vector<std::pair<uint32_t, uint32_t>> get_tf_both_directions_23mer_batch(const std::vector<std::string>& kmers) {
        std::vector<std::pair<uint32_t, uint32_t>> results;
        results.reserve(kmers.size());
        
        for (const auto& kmer : kmers) {
            results.push_back(get_tf_both_directions_23mer(kmer));
        }
        
        return results;
    }
    
    std::string get_reverse_complement_23mer(const std::string& kmer) {
        if (kmer.length() != 23) {
            return "";
        }
        
        std::string rev_kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
        uint64_t ukmer = get_dna23_bitset(kmer);
        uint64_t urev_kmer = reverseDNA(ukmer);
        get_bitset_dna23(urev_kmer, rev_kmer);
        
        return rev_kmer;
    }
    
    std::string get_23mer_statistics() {
        if (is_13mer_mode) {
            return "Not in 23-mer mode";
        }
        
        std::ostringstream stats;
        stats << "23-mer Index Statistics:\n";
        stats << "Total k-mers: " << n_kmers << "\n";
        stats << "Total reads: " << n_reads << "\n";
        stats << "AIndex loaded: " << (aindex_loaded ? "Yes" : "No") << "\n";
        stats << "Reads loaded: " << (reads != nullptr ? "Yes" : "No") << "\n";
        stats << "Hash map size: " << (hash_map ? hash_map->n : 0) << "\n";
        
        return stats.str();
    }
};

namespace py = pybind11;

PYBIND11_MODULE(aindex_cpp, m) {
    m.doc() = R"pbdoc(
        AIndex C++ Extension - K-mer Indexing and Querying Library
        
        High-performance C++ library for k-mer indexing and term frequency analysis.
        Supports both 13-mer and 23-mer indexing with memory-mapped file access
        and efficient hash-based lookups.
        
        Features:
        - Fast k-mer frequency counting and lookup
        - Memory-mapped file access for large datasets  
        - Support for both 13-mer and 23-mer k-mers
        - Reverse complement analysis
        - Read position tracking and retrieval
        - Batch processing capabilities
    )pbdoc";
    
    // Wrap the AindexWrapper class
    py::class_<AindexWrapper>(m, "AindexWrapper", R"pbdoc(
        Core AIndex wrapper class for k-mer indexing and querying.
        
        This class provides a high-level interface to the AIndex C++ library,
        supporting both 13-mer and 23-mer k-mer analysis with efficient
        memory-mapped file access and hash-based lookups.
    )pbdoc")
        .def(py::init<>(), R"pbdoc(
            Initialize a new AIndex wrapper instance.
            
            Creates an empty wrapper that can be configured to load
            either 13-mer or 23-mer indices using the load methods.
        )pbdoc")
        
        // =====================================================================
        // CATEGORY: Basic Index Loading Methods
        // =====================================================================
        .def("load", &AindexWrapper::load, R"pbdoc(
            [BASIC LOADING] Load complete 23-mer index from individual files.
            
            Loads a complete 23-mer index from the component files: hash file,
            term frequency file, k-mer binary file, and k-mer text file.
            
            Args:
                hash_file (str): Path to the hash index file (.pf)
                tf_file (str): Path to term frequency file (.tf.bin)  
                kmers_bin_file (str): Path to k-mer binary file (.kmers.bin)
                kmers_text_file (str): Path to k-mer text file (.kmers)
                
            Example:
                >>> wrapper = AindexWrapper()
                >>> wrapper.load("data.pf", "data.tf.bin", "data.kmers.bin", "data.kmers")
        )pbdoc")
        .def("load_hash_file", &AindexWrapper::load_hash_file, R"pbdoc(
            [BASIC LOADING] Load hash index and term frequency data.
            
            Loads the core hash index and term frequency data needed for
            k-mer frequency lookups. This is a lighter alternative to
            full index loading when position data is not needed.
            
            Args:
                hash_file (str): Path to the hash index file (.pf)
                tf_file (str): Path to term frequency file (.tf.bin)
                kmers_bin_file (str): Path to k-mer binary file (.kmers.bin) 
                kmers_text_file (str): Path to k-mer text file (.kmers)
                
            Example:
                >>> wrapper = AindexWrapper()
                >>> wrapper.load_hash_file("data.pf", "data.tf.bin", "data.kmers.bin", "data.kmers")
        )pbdoc")
        .def("load_reads", &AindexWrapper::load_reads, R"pbdoc(
            [BASIC LOADING] Load reads data from file.
            
            Loads read sequence data for position-based queries and read retrieval.
            The reads file contains the original sequences that were indexed.
            
            Args:
                reads_file (str): Path to the reads file
                
            Example:
                >>> wrapper.load_reads("sequences.fasta")
        )pbdoc")
        .def("load_reads_index", &AindexWrapper::load_reads_index, R"pbdoc(
            [BASIC LOADING] Load reads index for position tracking.
            
            Loads a pre-computed reads index that enables efficient mapping
            from k-mer positions back to specific reads and their coordinates.
            
            Args:
                reads_index_file (str): Path to the reads index file
                
            Example:
                >>> wrapper.load_reads_index("reads.index")
        )pbdoc")
        .def("load_reads_in_memory", &AindexWrapper::load_reads_in_memory, R"pbdoc(
            [BASIC LOADING] Load reads file completely into memory.
            
            Loads the entire reads file into memory for faster access.
            Use this for smaller datasets or when memory is not a constraint
            and maximum query speed is desired.
            
            Args:
                reads_file (str): Path to the reads file
                
            Example:
                >>> wrapper.load_reads_in_memory("small_dataset.fasta")
        )pbdoc")
        .def("load_aindex", &AindexWrapper::load_aindex, R"pbdoc(
            [BASIC LOADING] Load AIndex position data.
            
            Loads the AIndex position mapping data that enables efficient
            retrieval of k-mer positions within reads. Required for
            position-based queries and read context analysis.
            
            Args:
                aindex_file (str): Path to the AIndex file
                
            Example:
                >>> wrapper.load_aindex("data.aindex")
        )pbdoc")
        .def("load_13mer_index", &AindexWrapper::load_13mer_index, R"pbdoc(
            [13-MER LOADING] Load 13-mer hash index and term frequencies.
            
            Loads a 13-mer specific index consisting of hash mappings and
            term frequency data. 13-mer mode provides faster lookups for
            shorter k-mers with reduced memory usage.
            
            Args:
                hash_file (str): Path to 13-mer hash file
                tf_file (str): Path to 13-mer term frequency file
                
            Example:
                >>> wrapper.load_13mer_index("data_13mer.hash", "data_13mer.tf")
        )pbdoc")
        .def("load_13mer_aindex", &AindexWrapper::load_13mer_aindex, R"pbdoc(
            [13-MER LOADING] Load 13-mer position index from component files.
            
            Loads the 13-mer position index from its component files: position
            data, index mappings, and indices arrays. Enables position-based
            queries for 13-mer k-mers.
            
            Args:
                index_file (str): Path to index mappings file (.index.bin)
                indices_file (str): Path to indices array file (.indices.bin)
                
            Example:
                >>> wrapper.load_13mer_aindex("data.index.bin", "data.indices.bin")
        )pbdoc")
        
        // =====================================================================
        // CATEGORY: Convenient Prefix-Based Loading Methods  
        // =====================================================================
        .def("load_from_prefix_23mer", &AindexWrapper::load_from_prefix_23mer, R"pbdoc(
            [CONVENIENT LOADING] Load complete 23-mer index using file prefix.
            
            Automatically constructs file paths from a common prefix and loads
            a complete 23-mer index. This is the recommended method for loading
            23-mer indices as it handles all required files automatically.
            
            Args:
                prefix (str): Common file prefix (e.g., "data" for "data.pf", "data.tf.bin", etc.)
                reads_file (str, optional): Path to reads file to load alongside index
                
            File pattern:
                {prefix}.pf, {prefix}.tf.bin, {prefix}.kmers.bin, {prefix}.kmers
                
            Example:
                >>> wrapper = AindexWrapper()
                >>> wrapper.load_from_prefix_23mer("my_dataset", "reads.fasta")
        )pbdoc",
             py::arg("prefix"), py::arg("reads_file") = "")
        .def("load_aindex_from_prefix_23mer", &AindexWrapper::load_aindex_from_prefix_23mer, R"pbdoc(
            [CONVENIENT LOADING] Load 23-mer AIndex with position data using prefix.
            
            Loads a complete 23-mer index including position data for read context
            queries. Uses automatic file path construction and applies TF filtering
            for memory efficiency.
            
            Args:
                prefix (str): Common file prefix for index files
                max_tf (int): Maximum term frequency threshold for filtering
                reads_file (str, optional): Path to reads file to load
                
            File pattern:
                {prefix}.aindex, plus standard 23-mer files
                
            Example:
                >>> wrapper = AindexWrapper()
                >>> wrapper.load_aindex_from_prefix_23mer("dataset", max_tf=1000)
        )pbdoc",
             py::arg("prefix"), py::arg("max_tf"), py::arg("reads_file") = "")
        .def("load_from_prefix_13mer", &AindexWrapper::load_from_prefix_13mer, R"pbdoc(
            [CONVENIENT LOADING] Load complete 13-mer index using file prefix.
            
            Automatically constructs file paths and loads a complete 13-mer index.
            13-mer indices are more memory efficient and provide faster lookups
            for shorter k-mer analysis.
            
            Args:
                prefix (str): Common file prefix (e.g., "data13" for "data13.hash", etc.)
                reads_file (str, optional): Path to reads file to load
                
            File pattern:
                {prefix}.hash, {prefix}.tf.bin
                
            Example:
                >>> wrapper = AindexWrapper()
                >>> wrapper.load_from_prefix_13mer("my_13mer_data", "sequences.fasta")
        )pbdoc",
             py::arg("prefix"), py::arg("reads_file") = "")
        .def("load_aindex_from_prefix_13mer", &AindexWrapper::load_aindex_from_prefix_13mer, R"pbdoc(
            [CONVENIENT LOADING] Load 13-mer AIndex with position data using prefix.
            
            Loads a complete 13-mer index including position data for read context
            queries. Provides the full functionality of 13-mer analysis including
            position tracking and read retrieval.
            
            Args:
                prefix (str): Common file prefix for position files  
                reads_file (str, optional): Path to reads file to load
                
            File pattern:
                {prefix}.index.bin, {prefix}.indices.bin
                
            Example:
                >>> wrapper = AindexWrapper()
                >>> wrapper.load_aindex_from_prefix_13mer("pos_data", "reads.fasta")
        )pbdoc",
             py::arg("prefix"), py::arg("reads_file") = "")
        
        // =====================================================================
        // CATEGORY: Term Frequency and Hash Queries
        // =====================================================================
        .def("get_tf_values", &AindexWrapper::get_tf_values, R"pbdoc(
            [TF QUERIES] Get term frequency values for multiple k-mers.
            
            Retrieves term frequency (TF) values for a batch of k-mers.
            Term frequency represents how many times each k-mer appears
            in the indexed dataset.
            
            Args:
                kmers (List[str]): List of k-mer sequences to query
                
            Returns:
                List[int]: Term frequency values for each k-mer
                
            Example:
                >>> tf_values = wrapper.get_tf_values(["ATCGATCGATCG", "GCTAGCTAGCTA"])
                >>> print(f"TF values: {tf_values}")
        )pbdoc")
        .def("get_tf_value", &AindexWrapper::get_tf_value, R"pbdoc(
            [TF QUERIES] Get term frequency value for a single k-mer.
            
            Retrieves the term frequency (TF) value for a single k-mer.
            Returns 0 if the k-mer is not found in the index.
            
            Args:
                kmer (str): K-mer sequence to query
                
            Returns:
                int: Term frequency value for the k-mer
                
            Example:
                >>> tf = wrapper.get_tf_value("ATCGATCGATCG")
                >>> print(f"K-mer frequency: {tf}")
        )pbdoc")
        .def("get_hash_values", &AindexWrapper::get_hash_values, R"pbdoc(
            [HASH QUERIES] Get hash values for multiple k-mers.
            
            Retrieves hash index values for a batch of k-mers. Hash values
            are internal identifiers used for efficient k-mer lookup.
            
            Args:
                kmers (List[str]): List of k-mer sequences to hash
                
            Returns:
                List[int]: Hash values for each k-mer
                
            Example:
                >>> hashes = wrapper.get_hash_values(["ATCGATCGATCG", "GCTAGCTAGCTA"])
                >>> print(f"Hash values: {hashes}")
        )pbdoc")
        .def("get_hash_value", &AindexWrapper::get_hash_value, R"pbdoc(
            [HASH QUERIES] Get hash value for a single k-mer.
            
            Retrieves the hash index value for a single k-mer.
            Hash values are used internally for efficient k-mer indexing.
            
            Args:
                kmer (str): K-mer sequence to hash
                
            Returns:
                int: Hash value for the k-mer
                
            Example:
                >>> hash_val = wrapper.get_hash_value("ATCGATCGATCG")
                >>> print(f"Hash: {hash_val}")
        )pbdoc")
        .def("get_kid_by_kmer", &AindexWrapper::get_kid_by_kmer, R"pbdoc(
            [HASH QUERIES] Get k-mer ID by k-mer sequence.
            
            Retrieves the internal k-mer ID (kid) for a given k-mer sequence.
            The k-mer ID is used internally for efficient k-mer indexing and lookup.
            
            Args:
                kmer (str): K-mer sequence to get ID for
                
            Returns:
                int: K-mer ID for the k-mer
                
            Example:
                >>> kid = wrapper.get_kid_by_kmer("ATCGATCGATCG")
                >>> print(f"K-mer ID: {kid}")
        )pbdoc")
        .def("get_kmer_by_kid", &AindexWrapper::get_kmer_by_kid, R"pbdoc(
            [HASH QUERIES] Get k-mer sequence by k-mer ID.
            
            Retrieves the k-mer sequence for a given internal k-mer ID (kid).
            This is the reverse operation of get_kid_by_kmer.
            
            Args:
                kid (int): K-mer ID to get sequence for
                
            Returns:
                str: K-mer sequence for the ID
                
            Example:
                >>> kmer = wrapper.get_kmer_by_kid(42)
                >>> print(f"K-mer for ID 42: {kmer}")
        )pbdoc")
        .def("get_strand", &AindexWrapper::get_strand, R"pbdoc(
            [HASH QUERIES] Get strand information for a k-mer.
            
            Determines the strand (forward, reverse, or not found) for a given k-mer
            by checking if it exists in forward or reverse complement form.
            
            Args:
                kmer (str): K-mer sequence to check strand for
                
            Returns:
                int: Strand value (0=not found, 1=forward, 2=reverse)
                
            Example:
                >>> strand = wrapper.get_strand("ATCGATCGATCG")
                >>> print(f"Strand: {strand}")
        )pbdoc")
        .def("get_kmer_info", &AindexWrapper::get_kmer_info, R"pbdoc(
            [HASH QUERIES] Get comprehensive k-mer information by k-mer ID.
            
            Retrieves detailed information about a k-mer including its term frequency,
            forward sequence, and reverse complement sequence using its k-mer ID.
            
            Args:
                kid (int): K-mer ID to get information for
                
            Returns:
                Tuple[int, str, str]: (term_frequency, kmer, reverse_complement_kmer)
                
            Example:
                >>> tf, kmer, rkmer = wrapper.get_kmer_info(42)
                >>> print(f"TF: {tf}, K-mer: {kmer}, Rev-comp: {rkmer}")
        )pbdoc")
        .def("get_rid", &AindexWrapper::get_rid, R"pbdoc(
            [POSITION QUERIES] Get read ID by position.
            
            Retrieves the read ID that contains the specified position in the
            concatenated reads data. Requires position index to be loaded.
            
            Args:
                pos (int): Position to query
                
            Returns:
                int: Read ID containing the position
                
            Example:
                >>> rid = wrapper.get_rid(1000)
                >>> print(f"Position 1000 is in read: {rid}")
        )pbdoc")
        .def("get_start", &AindexWrapper::get_start, R"pbdoc(
            [POSITION QUERIES] Get start position of read containing given position.
            
            Retrieves the start position of the read that contains the specified
            position in the concatenated reads data.
            
            Args:
                pos (int): Position to query
                
            Returns:
                int: Start position of the read containing pos
                
            Example:
                >>> start = wrapper.get_start(1000)
                >>> print(f"Read containing position 1000 starts at: {start}")
        )pbdoc")
        
        // =====================================================================
        // CATEGORY: Read and Position Queries
        // =====================================================================
        .def("get_read_by_rid", &AindexWrapper::get_read_by_rid, R"pbdoc(
            [READ QUERIES] Get read sequence by read ID.
            
            Retrieves the complete sequence of a read using its unique read ID.
            Requires reads data to be loaded via load_reads() or similar method.
            
            Args:
                rid (int): Read ID to retrieve
                
            Returns:
                str: Read sequence
                
            Example:
                >>> read_seq = wrapper.get_read_by_rid(42)
                >>> print(f"Read 42: {read_seq}")
        )pbdoc")
        .def("get_read", &AindexWrapper::get_read, R"pbdoc(
            [READ QUERIES] Get read sequence by position range.
            
            Retrieves a read sequence from the specified start and end positions
            within the concatenated reads data. Optionally returns the reverse
            complement if requested.
            
            Args:
                start (int): Start position in reads data
                end (int): End position in reads data  
                revcomp (bool): Return reverse complement if True
                
            Returns:
                str: Read sequence (or reverse complement)
                
            Example:
                >>> seq = wrapper.get_read(1000, 1023, revcomp=False)
                >>> rev_seq = wrapper.get_read(1000, 1023, revcomp=True)
        )pbdoc", py::arg("start"), py::arg("end"), py::arg("revcomp") = false)
        .def("get_reads_se_by_kmer", &AindexWrapper::get_reads_se_by_kmer, R"pbdoc(
            [READ QUERIES] Get reads containing a specific k-mer.
            
            Finds all reads that contain the specified k-mer and returns
            information about their positions and context. Useful for
            analyzing k-mer distribution across reads.
            
            Args:
                kmer (str): K-mer sequence to search for
                
            Returns:
                List: Information about reads containing the k-mer
                
            Example:
                >>> reads_info = wrapper.get_reads_se_by_kmer("ATCGATCGATCG")
                >>> print(f"Found in {len(reads_info)} reads")
        )pbdoc")
        .def("get_positions", 
             static_cast<std::vector<uint64_t>(AindexWrapper::*)(const std::string&)>(&AindexWrapper::get_positions), 
             R"pbdoc(
            [POSITION QUERIES] Get all positions where a k-mer appears.
            
            Retrieves all positions in the indexed data where the specified
            k-mer occurs. Requires AIndex position data to be loaded.
            
            Args:
                kmer (str): K-mer sequence to locate
                
            Returns:
                List[int]: List of positions where k-mer appears
                
            Example:
                >>> positions = wrapper.get_positions("ATCGATCGATCG")
                >>> print(f"K-mer found at positions: {positions}")
        )pbdoc")
        
        // =====================================================================
        // CATEGORY: Index Statistics and Metadata
        // =====================================================================
        .def("get_hash_size", &AindexWrapper::get_hash_size, R"pbdoc(
            [METADATA] Get size of the hash index.
            
            Returns the total size of the hash index structure used
            for k-mer lookups. Useful for memory usage analysis.
            
            Returns:
                int: Size of hash index
                
            Example:
                >>> size = wrapper.get_hash_size()
                >>> print(f"Hash index size: {size}")
        )pbdoc")
        .def("get_reads_size", &AindexWrapper::get_reads_size, R"pbdoc(
            [METADATA] Get total size of reads data.
            
            Returns the total size of the loaded reads data in characters/bases.
            Useful for understanding dataset scale and memory usage.
            
            Returns:
                int: Total size of reads data
                
            Example:
                >>> size = wrapper.get_reads_size()
                >>> print(f"Total reads size: {size} bases")
        )pbdoc")
        
        // =====================================================================
        // CATEGORY: Instance Properties (Read/Write)
        // =====================================================================
        .def_readwrite("aindex_loaded", &AindexWrapper::aindex_loaded, R"pbdoc(
            [PROPERTY] Boolean flag indicating if AIndex position data is loaded.
            
            True if position index data has been successfully loaded,
            enabling position-based queries and read context analysis.
        )pbdoc")
        .def_readwrite("n_reads", &AindexWrapper::n_reads, R"pbdoc(
            [PROPERTY] Number of reads in the loaded dataset.
            
            Total count of individual reads/sequences in the indexed dataset.
            Updated when reads data is loaded.
        )pbdoc")
        .def_readwrite("n_kmers", &AindexWrapper::n_kmers, R"pbdoc(
            [PROPERTY] Number of unique k-mers in the index.
            
            Total count of unique k-mers that have been indexed.
            This represents the vocabulary size of the k-mer index.
        )pbdoc")
        .def_readwrite("reads_size", &AindexWrapper::reads_size, R"pbdoc(
            [PROPERTY] Total size of reads data in characters.
            
            Total number of characters/bases across all loaded reads.
            Useful for memory usage calculations and dataset statistics.
        )pbdoc")
        
        // =====================================================================
        // CATEGORY: Debugging and Development Tools
        // =====================================================================
        .def("debug_kmer_tf_values", &AindexWrapper::debug_kmer_tf_values, R"pbdoc(
            [DEBUG] Debug k-mer term frequency values.
            
            Diagnostic function for debugging k-mer frequency calculations.
            Provides detailed information about internal TF value computation.
            
            Example:
                >>> wrapper.debug_kmer_tf_values()
        )pbdoc")
        .def("get_index_info", &AindexWrapper::get_index_info, R"pbdoc(
            [DEBUG] Get comprehensive index statistics and information.
            
            Returns detailed information about the loaded index including
            memory usage, k-mer counts, file sizes, and loading status.
            
            Returns:
                str: Detailed index information and statistics
                
            Example:
                >>> info = wrapper.get_index_info()
                >>> print(info)
        )pbdoc")
        
        // =====================================================================
        // CATEGORY: Specialized 13-mer Analysis Methods
        // =====================================================================
        .def("get_total_tf_value_13mer", &AindexWrapper::get_total_tf_value_13mer, R"pbdoc(
            [13-MER ANALYSIS] Get combined TF value for 13-mer and its reverse complement.
            
            Calculates the total term frequency by summing the TF values of
            the forward k-mer and its reverse complement. Useful for strand-
            agnostic k-mer frequency analysis.
            
            Args:
                kmer (str): 13-mer sequence (must be length 13)
                
            Returns:
                int: Combined TF value (forward + reverse complement)
                
            Example:
                >>> total_tf = wrapper.get_total_tf_value_13mer("ATCGATCGATCGA")
                >>> print(f"Total 13-mer frequency: {total_tf}")
        )pbdoc")
        .def("get_total_tf_values_13mer", &AindexWrapper::get_total_tf_values_13mer, R"pbdoc(
            [13-MER ANALYSIS] Get combined TF values for multiple 13-mers and their reverse complements.
            
            Batch version of get_total_tf_value_13mer(). Efficiently calculates
            combined forward and reverse complement frequencies for multiple 13-mers.
            
            Args:
                kmers (List[str]): List of 13-mer sequences
                
            Returns:
                List[int]: Combined TF values for each 13-mer
                
            Example:
                >>> total_tfs = wrapper.get_total_tf_values_13mer(["ATCGATCGATCGA", "GCTAGCTAGCTAG"])
                >>> print(f"Total frequencies: {total_tfs}")
        )pbdoc")
        .def("get_tf_both_directions_13mer", &AindexWrapper::get_tf_both_directions_13mer, R"pbdoc(
            [13-MER ANALYSIS] Get separate TF values for 13-mer forward and reverse complement.
            
            Returns the term frequency values for both the forward k-mer and
            its reverse complement as separate values. Useful for strand-specific
            analysis and understanding k-mer directionality.
            
            Args:
                kmer (str): 13-mer sequence (must be length 13)
                
            Returns:
                Tuple[int, int]: (forward_tf, reverse_complement_tf)
                
            Example:
                >>> fw_tf, rv_tf = wrapper.get_tf_both_directions_13mer("ATCGATCGATCGA")
                >>> print(f"Forward: {fw_tf}, Reverse: {rv_tf}")
        )pbdoc")
        .def("get_tf_both_directions_13mer_batch", &AindexWrapper::get_tf_both_directions_13mer_batch, R"pbdoc(
            [13-MER ANALYSIS] Get directional TF values for multiple 13-mers.
            
            Batch version of get_tf_both_directions_13mer(). Efficiently retrieves
            forward and reverse complement TF values for multiple 13-mers.
            
            Args:
                kmers (List[str]): List of 13-mer sequences
                
            Returns:
                List[Tuple[int, int]]: List of (forward_tf, reverse_tf) pairs
                
            Example:
                >>> results = wrapper.get_tf_both_directions_13mer_batch(["ATCGATCGATCGA", "GCTAGCTAGCTAG"])
                >>> for fw, rv in results:
                >>>     print(f"Forward: {fw}, Reverse: {rv}")
        )pbdoc")
        .def("get_reverse_complement_13mer", &AindexWrapper::get_reverse_complement_13mer, R"pbdoc(
            [13-MER ANALYSIS] Get reverse complement of a 13-mer sequence.
            
            Computes the reverse complement of a 13-mer k-mer using optimized
            bit operations. Part of the 13-mer analysis toolkit.
            
            Args:
                kmer (str): 13-mer sequence (must be length 13)
                
            Returns:
                str: Reverse complement sequence
                
            Example:
                >>> rev_comp = wrapper.get_reverse_complement_13mer("ATCGATCGATCGA")
                >>> print(f"Original: ATCGATCGATCGA, Rev comp: {rev_comp}")
        )pbdoc")
        .def("get_13mer_statistics", &AindexWrapper::get_13mer_statistics, R"pbdoc(
            [13-MER ANALYSIS] Get comprehensive statistics about the 13-mer index.
            
            Returns detailed statistics about the loaded 13-mer index including
            total k-mer count, read count, memory usage, and loading status.
            
            Returns:
                str: Formatted statistics string
                
            Example:
                >>> stats = wrapper.get_13mer_statistics()
                >>> print(stats)
        )pbdoc")
        .def("get_13mer_tf_array", &AindexWrapper::get_13mer_tf_array, R"pbdoc(
            [13-MER ANALYSIS] Get direct access to the 13-mer term frequency array.
            
            Provides low-level access to the internal TF array for 13-mers.
            Advanced users can use this for custom analysis and optimization.
            
            Returns:
                Array: Direct reference to 13-mer TF array
                
            Example:
                >>> tf_array = wrapper.get_13mer_tf_array()
                >>> print(f"Array size: {len(tf_array)}")
        )pbdoc")
        .def("get_tf_by_index_13mer", &AindexWrapper::get_tf_by_index_13mer, R"pbdoc(
            [13-MER ANALYSIS] Get TF value by direct array index for 13-mers.
            
            Retrieves term frequency value using direct array indexing.
            This is the fastest method for TF lookup when you know the
            internal index of the k-mer.
            
            Args:
                index (int): Internal array index
                
            Returns:
                int: Term frequency value at index
                
            Example:
                >>> tf = wrapper.get_tf_by_index_13mer(1234)
                >>> print(f"TF at index 1234: {tf}")
        )pbdoc")
        .def("get_positions_13mer", 
             static_cast<std::vector<uint64_t>(AindexWrapper::*)(const std::string&)>(&AindexWrapper::get_positions_13mer), 
             R"pbdoc(
            [13-MER ANALYSIS] Get positions for 13-mers using the position index.
            
            Retrieves all positions where 13-mer k-mers appear in the indexed
            data using the specialized 13-mer position index for optimal performance.
            
            Args:
                kmer (str): 13-mer sequence to locate
                
            Returns:
                List[int]: Positions where 13-mer appears
                
            Example:
                >>> positions = wrapper.get_positions_13mer("ATCGATCGATCGA")
                >>> print(f"13-mer found at: {positions}")
        )pbdoc")
        
        // =====================================================================
        // CATEGORY: Specialized 23-mer Analysis Methods
        // =====================================================================
        .def("get_tf_values_23mer", &AindexWrapper::get_tf_values_23mer, R"pbdoc(
            [23-MER ANALYSIS] Get term frequency values for 23-mer k-mers.
            
            Optimized method for retrieving term frequency values specifically
            for 23-mer k-mers. Uses the 23-mer index for efficient lookup.
            
            Args:
                kmers (List[str]): List of 23-mer sequences (must be length 23)
                
            Returns:
                List[int]: Term frequency values for each 23-mer
                
            Example:
                >>> tf_vals = wrapper.get_tf_values_23mer(["ATCGATCGATCGATCGATCGATC", "GCTAGCTAGCTAGCTAGCTAGCT"])
                >>> print(f"23-mer TF values: {tf_vals}")
        )pbdoc")
        .def("get_total_tf_value_23mer", &AindexWrapper::get_total_tf_value_23mer, R"pbdoc(
            [23-MER ANALYSIS] Get combined TF value for 23-mer and its reverse complement.
            
            Calculates the total term frequency by summing the TF values of
            the forward k-mer and its reverse complement. Useful for strand-
            agnostic k-mer frequency analysis in longer sequences.
            
            Args:
                kmer (str): 23-mer sequence (must be length 23)
                
            Returns:
                int: Combined TF value (forward + reverse complement)
                
            Example:
                >>> total_tf = wrapper.get_total_tf_value_23mer("ATCGATCGATCGATCGATCGATC")
                >>> print(f"Total 23-mer frequency: {total_tf}")
        )pbdoc")
        .def("get_total_tf_values_23mer", &AindexWrapper::get_total_tf_values_23mer, R"pbdoc(
            [23-MER ANALYSIS] Get combined TF values for multiple 23-mers and their reverse complements.
            
            Batch version of get_total_tf_value_23mer(). Efficiently calculates
            combined forward and reverse complement frequencies for multiple 23-mers.
            
            Args:
                kmers (List[str]): List of 23-mer sequences
                
            Returns:
                List[int]: Combined TF values for each 23-mer
                
            Example:
                >>> total_tfs = wrapper.get_total_tf_values_23mer(["ATCGATCGATCGATCGATCGATC", "GCTAGCTAGCTAGCTAGCTAGCT"])
                >>> print(f"Total frequencies: {total_tfs}")
        )pbdoc")
        .def("get_tf_both_directions_23mer", &AindexWrapper::get_tf_both_directions_23mer, R"pbdoc(
            [23-MER ANALYSIS] Get separate TF values for 23-mer forward and reverse complement.
            
            Returns the term frequency values for both the forward k-mer and
            its reverse complement as separate values. Useful for strand-specific
            analysis and understanding k-mer directionality in longer sequences.
            
            Args:
                kmer (str): 23-mer sequence (must be length 23)
                
            Returns:
                Tuple[int, int]: (forward_tf, reverse_complement_tf)
                
            Example:
                >>> fw_tf, rv_tf = wrapper.get_tf_both_directions_23mer("ATCGATCGATCGATCGATCGATC")
                >>> print(f"Forward: {fw_tf}, Reverse: {rv_tf}")
        )pbdoc")
        .def("get_tf_both_directions_23mer_batch", &AindexWrapper::get_tf_both_directions_23mer_batch, R"pbdoc(
            [23-MER ANALYSIS] Get directional TF values for multiple 23-mers.
            
            Batch version of get_tf_both_directions_23mer(). Efficiently retrieves
            forward and reverse complement TF values for multiple 23-mers.
            
            Args:
                kmers (List[str]): List of 23-mer sequences
                
            Returns:
                List[Tuple[int, int]]: List of (forward_tf, reverse_tf) pairs
                
            Example:
                >>> results = wrapper.get_tf_both_directions_23mer_batch(["ATCGATCGATCGATCGATCGATC", "GCTAGCTAGCTAGCTAGCTAGCT"])
                >>> for fw, rv in results:
                >>>     print(f"Forward: {fw}, Reverse: {rv}")
        )pbdoc")
        .def("get_reverse_complement_23mer", &AindexWrapper::get_reverse_complement_23mer, R"pbdoc(
            [23-MER ANALYSIS] Get reverse complement of a 23-mer sequence.
            
            Computes the reverse complement of a 23-mer k-mer using optimized
            bit operations. Part of the 23-mer analysis toolkit.
            
            Args:
                kmer (str): 23-mer sequence (must be length 23)
                
            Returns:
                str: Reverse complement sequence
                
            Example:
                >>> rev_comp = wrapper.get_reverse_complement_23mer("ATCGATCGATCGATCGATCGATC")
                >>> print(f"Original: ATCGATCGATCGATCGATCGATC")
                >>> print(f"Rev comp: {rev_comp}")
        )pbdoc")
        .def("get_23mer_statistics", &AindexWrapper::get_23mer_statistics, R"pbdoc(
            [23-MER ANALYSIS] Get comprehensive statistics about the 23-mer index.
            
            Returns detailed statistics about the loaded 23-mer index including
            total k-mer count, read count, memory usage, and loading status.
            
            Returns:
                str: Formatted statistics string
                
            Example:
                >>> stats = wrapper.get_23mer_statistics()
                >>> print(stats)
        )pbdoc");
}
