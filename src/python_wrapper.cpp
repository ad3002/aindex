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
    uint32_t* tf_array_13mer = nullptr;  // Use uint32_t to match file format
    static const uint32_t TOTAL_13MERS = 67108864; // 4^13
    
    // 13-mer AIndex support (positions)
    uint64_t *positions_13mer = nullptr;
    uint64_t *indices_13mer = nullptr;
    uint64_t n_13mer = 0;
    uint32_t max_tf_13mer = 0;
    uint64_t indices_length_13mer = 0;
    
public:
    bool aindex_loaded = false;
    PHASH_MAP *hash_map;
    uint64_t n_reads = 0;
    uint64_t n_kmers = 0;
    
    uint64_t reads_size = 0;
    char *reads = nullptr;

    std::unordered_map<uint64_t, uint32_t> start2rid;
    std::unordered_map<uint64_t, uint64_t> start2end;
    std::vector<uint64_t> start_positions;

    IntervalTree pos_intervalTree;
    
    AindexWrapper() {}

    ~AindexWrapper() {
        if (positions != nullptr) munmap(positions, n*sizeof(uint64_t));
        if (indices != nullptr) munmap(indices, indices_length);
        if (reads != nullptr) munmap(reads, reads_size);
        if (tf_array_13mer != nullptr) munmap(tf_array_13mer, TOTAL_13MERS * sizeof(uint32_t));
        if (positions_13mer != nullptr) munmap(positions_13mer, indices_length_13mer);
        if (indices_13mer != nullptr) munmap(indices_13mer, (TOTAL_13MERS + 1) * sizeof(uint64_t));
        if (positions_13mer != nullptr) munmap(positions_13mer, n_13mer*sizeof(uint64_t));
        if (indices_13mer != nullptr) munmap(indices_13mer, indices_length_13mer);

        delete hash_map;

        reads = nullptr;
        indices = nullptr;
        positions = nullptr;
        tf_array_13mer = nullptr;
        positions_13mer = nullptr;
        indices_13mer = nullptr;
    }

    void load(std::string hash_filename, std::string tf_file, std::string kmers_bin_filename, std::string kmers_text_filename){
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
        // Memory map reads
        emphf::logger() << "Memory mapping reads file..." << std::endl;
        std::ifstream fout(reads_file, std::ios::in | std::ios::binary);
        fout.seekg(0, std::ios::end);
        uint64_t length = fout.tellg();
        fout.close();

        FILE* in = std::fopen(reads_file.c_str(), "rb");
        reads = (char*)mmap(NULL, length, PROT_READ|PROT_WRITE, MAP_PRIVATE, fileno(in), 0);
        if (reads == nullptr) {
            std::cerr << "Failed position loading" << std::endl;
            exit(10);
        }
        fclose(in);

        reads_size = length;

        emphf::logger() << "\tbuilding start pos index over reads: " << std::endl;
        std::string index_file = reads_file.substr(0, reads_file.find_last_of(".")) + ".ridx";
        load_reads_index(index_file);
        emphf::logger() << "\tDone" << std::endl;
    }

    void load_reads_in_memory(std::string reads_file) {
        // Load reads into memory
        emphf::logger() << "Loading reads file into memory..." << std::endl;
        std::ifstream fin(reads_file, std::ios::in | std::ios::binary);
        if (!fin) {
            std::cerr << "Failed to open file" << std::endl;
            exit(1);
        }

        fin.seekg(0, std::ios::end);
        uint64_t length = fin.tellg();
        fin.seekg(0, std::ios::beg);

        reads = new char[length];
        fin.read(reads, length);
        fin.close();

        reads_size = length;

        emphf::logger() << "\tbuilding start pos index over reads: " << std::endl;
        std::string index_file = reads_file.substr(0, reads_file.find_last_of(".")) + ".ridx";
        load_reads_index(index_file);
        emphf::logger() << "\tDone" << std::endl;
    }

    void load_aindex(std::string pos_file, std::string index_file, std::string indices_file, uint32_t _max_tf) {
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
        
        tf_array_13mer = (uint32_t*)mmap(NULL, TOTAL_13MERS * sizeof(uint32_t), 
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
    
    void load_13mer_aindex(const std::string& pos_file, const std::string& index_file, const std::string& indices_file) {
        emphf::logger() << "Loading 13-mer AIndex files..." << std::endl;
        emphf::logger() << "Pos file: " << pos_file << std::endl;
        emphf::logger() << "Index file: " << index_file << std::endl;
        emphf::logger() << "Indices file: " << indices_file << std::endl;
        
        // Load pos.bin file
        std::ifstream pos_in(pos_file, std::ios::in | std::ios::binary);
        if (!pos_in) {
            std::cerr << "Failed to open pos file: " << pos_file << std::endl;
            std::terminate();
        }
        pos_in.seekg(0, std::ios::end);
        uint64_t pos_length = pos_in.tellg();
        pos_in.close();
        
        n_13mer = pos_length / sizeof(uint64_t);
        emphf::logger() << "\tPositions length: " << pos_length << " (" << n_13mer << " positions)" << std::endl;
        
        FILE* pos_fp = std::fopen(pos_file.c_str(), "rb");
        if (!pos_fp) {
            std::cerr << "Failed to open pos file for mmap: " << pos_file << std::endl;
            std::terminate();
        }
        
        positions_13mer = (uint64_t*)mmap(NULL, pos_length, PROT_READ, MAP_SHARED, fileno(pos_fp), 0);
        if (positions_13mer == MAP_FAILED) {
            std::cerr << "Failed to mmap pos file" << std::endl;
            std::terminate();
        }
        fclose(pos_fp);
        
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
        
        emphf::logger() << "13-mer AIndex loaded successfully" << std::endl;
    }
    
    bool is_13mer(const std::string& kmer) const {
        return kmer.length() == 13;
    }
    
    bool is_23mer(const std::string& kmer) const {
        return kmer.length() == 23;
    }
    
    uint32_t get_tf_value_13mer(const std::string& kmer) {
        if (!is_13mer_mode || !is_13mer(kmer)) {
            return 0;
        }
        
        // Validate k-mer (only ATGC)
        for (char c : kmer) {
            if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
                return 0;
            }
        }
        
        emphf::stl_string_adaptor str_adapter;
        
        // Try forward direction
        uint64_t hash_id = hasher_13mer.lookup(kmer, str_adapter);
        if (hash_id < TOTAL_13MERS) {
            return tf_array_13mer[hash_id];
        }
        
        // Try reverse complement
        std::string rev_kmer = get_reverse_complement_13mer(kmer);
        hash_id = hasher_13mer.lookup(rev_kmer, str_adapter);
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

    std::vector<uint32_t> get_tf_values(std::vector<std::string> kmers) {
        std::vector<uint32_t> tf_values;
        
        // Auto-detect mode from first k-mer
        if (!kmers.empty()) {
            if (is_13mer(kmers[0]) && is_13mer_mode) {
                // Use 13-mer batch processing
                tf_values.reserve(kmers.size());
                
                emphf::stl_string_adaptor str_adapter;
                
                for (const auto& kmer : kmers) {
                    if (!is_13mer(kmer)) {
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
                    
                    // Try forward direction
                    uint64_t hash_id = hasher_13mer.lookup(kmer, str_adapter);
                    if (hash_id < TOTAL_13MERS) {
                        tf_values.push_back(tf_array_13mer[hash_id]);
                        continue;
                    }
                    
                    // Try reverse complement
                    std::string rev_kmer = get_reverse_complement_13mer(kmer);
                    hash_id = hasher_13mer.lookup(rev_kmer, str_adapter);
                    if (hash_id < TOTAL_13MERS) {
                        tf_values.push_back(tf_array_13mer[hash_id]);
                    } else {
                        tf_values.push_back(0);
                    }
                }
                return tf_values;
            } else if (is_23mer(kmers[0]) && hash_map != nullptr) {
                // Original 23-mer logic
                for (const auto& kmer : kmers) {
                    uint64_t kmer_id = hash_map->hasher.lookup(kmer, str_adapter);
                    tf_values.push_back(hash_map->tf_values[kmer_id]);
                }
                return tf_values;
            }
        }
        
        // Default: return zeros for unknown k-mer sizes
        tf_values.resize(kmers.size(), 0);
        return tf_values;
    }

    uint32_t get_tf_value(std::string kmer) {
        // Auto-detect mode based on k-mer length
        if (is_13mer(kmer) && is_13mer_mode) {
            return get_tf_value_13mer(kmer);
        } else if (is_23mer(kmer) && hash_map != nullptr) {
            return get_tf_value_23mer(kmer);
        }
        
        return 0;
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

    std::string get_read_by_rid(uint64_t rid) {
        if (start_positions.size() <= rid) {
            return "";
        }
        uint64_t start = start_positions[rid];
        uint64_t end = start2end[start];
        
        std::string read(reads + start, end - start + 1);
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

    uint64_t get_kmer_info(uint64_t kid, std::string& kmer, std::string& rkmer) {
        if (kid >= hash_map->n) {
            return 0;
        }
        uint64_t ukmer = hash_map->checker[kid];
        kmer = get_bitset_dna23(ukmer);
        uint64_t urev_kmer = reverseDNA(ukmer);
        rkmer = get_bitset_dna23(urev_kmer);
        return hash_map->tf_values[kid];
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

    std::vector<uint64_t> get_positions(const std::string& kmer) {
        // Route to appropriate function based on k-mer length
        if (is_13mer(kmer)) {
            return get_positions_13mer(kmer);
        } else if (is_23mer(kmer)) {
            // Use the original algorithm for 23-mers
            std::vector<uint64_t> r;
            if (!hash_map || !positions || !indices) {
                return r;
            }
            auto h1 = hash_map->get_pfid(kmer);
            for (uint64_t i=indices[h1]; i < indices[h1+1] && h1+1 < indices_length; ++i) {
                if (positions[i] == 0) {
                    continue;
                }
                r.push_back(positions[i]-1);
            }
            return r;
        } else {
            // Unsupported k-mer length
            return std::vector<uint64_t>();
        }
    }

    uint64_t get_hash_size() {
        return hash_map ? hash_map->n : 0;
    }

    uint64_t get_reads_size() {
        return n_reads;
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
            if (!is_13mer(kmer)) {
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
            
            // Try forward direction
            uint64_t hash_id = hasher_13mer.lookup(kmer, str_adapter);
            if (hash_id < TOTAL_13MERS) {
                tf_values.push_back(tf_array_13mer[hash_id]);
                continue;
            }
            
            // Try reverse complement
            std::string rev_kmer = get_reverse_complement_13mer(kmer);
            hash_id = hasher_13mer.lookup(rev_kmer, str_adapter);
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
                    total_count += tf_array_13mer[i];  // uint32_t automatically promotes to uint64_t
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
    
    std::vector<uint64_t> get_positions_13mer(const std::string& kmer) {
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
        
        // Try forward direction
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
            return result;
        }
        
        // Try reverse complement
        std::string rev_kmer = get_reverse_complement_13mer(kmer);
        hash_id = hasher_13mer.lookup(rev_kmer, str_adapter);
        if (hash_id < TOTAL_13MERS) {
            // Get positions for reverse complement
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
};

namespace py = pybind11;

PYBIND11_MODULE(aindex_cpp, m) {
    m.doc() = "Aindex C++ bindings using pybind11";
    
    // Wrap the AindexWrapper class
    py::class_<AindexWrapper>(m, "AindexWrapper")
        .def(py::init<>())
        
        // Loading functions with new signature
        .def("load", &AindexWrapper::load,
             "Load index from hash file, tf file, kmers bin file, and kmers text file")
        .def("load_hash_file", &AindexWrapper::load_hash_file,
             "Load hash file, tf file, kmers bin file, and kmers text file") 
        .def("load_reads", &AindexWrapper::load_reads,
             "Load reads file")
        .def("load_reads_index", &AindexWrapper::load_reads_index,
             "Load reads index file")
        .def("load_reads_in_memory", &AindexWrapper::load_reads_in_memory,
             "Load reads file into memory")
        .def("load_aindex", &AindexWrapper::load_aindex,
             "Load aindex file")
        .def("load_13mer_index", &AindexWrapper::load_13mer_index,
             "Load 13-mer index from hash file and tf file")
        .def("load_13mer_aindex", &AindexWrapper::load_13mer_aindex,
             "Load 13-mer position index from pos, index, and indices files")
        .def("load_13mer_aindex", &AindexWrapper::load_13mer_aindex,
             "Load 13-mer AIndex files")
        
        // Query functions
        .def("get_tf_values", &AindexWrapper::get_tf_values,
             "Get term frequency values for kmers")
        .def("get_tf_value", &AindexWrapper::get_tf_value,
             "Get term frequency value for a kmer")
        .def("get_hash_values", &AindexWrapper::get_hash_values,
             "Get hash values for kmers")
        .def("get_hash_value", &AindexWrapper::get_hash_value,
             "Get hash value for a kmer")
        .def("get_read_by_rid", &AindexWrapper::get_read_by_rid,
             "Get read by read ID")
        .def("get_read", &AindexWrapper::get_read,
             "Get read by start and end positions", py::arg("start"), py::arg("end"), py::arg("revcomp") = false)
        .def("get_reads_se_by_kmer", &AindexWrapper::get_reads_se_by_kmer,
             "Get reads containing a specific kmer")
        .def("get_kid_by_kmer", &AindexWrapper::get_kid_by_kmer,
             "Get kmer ID by kmer")
        .def("get_kmer_by_kid", &AindexWrapper::get_kmer_by_kid,
             "Get kmer by kmer ID")
        .def("get_strand", &AindexWrapper::get_strand,
             "Get strand for kmer")
        .def("get_kmer_info", &AindexWrapper::get_kmer_info,
             "Get kmer info by kmer ID")
        .def("get_rid", &AindexWrapper::get_rid,
             "Get read ID by position")
        .def("get_start", &AindexWrapper::get_start,
             "Get start position by position")
        .def("get_positions", &AindexWrapper::get_positions,
             "Get positions for kmer")
        .def("get_hash_size", &AindexWrapper::get_hash_size,
             "Get hash size")
        .def("get_reads_size", &AindexWrapper::get_reads_size,
             "Get reads size")
        
        // Properties
        .def_readwrite("aindex_loaded", &AindexWrapper::aindex_loaded)
        .def_readwrite("n_reads", &AindexWrapper::n_reads)
        .def_readwrite("n_kmers", &AindexWrapper::n_kmers)
        .def_readwrite("reads_size", &AindexWrapper::reads_size)
        
        // Debug function
        .def("debug_kmer_tf_values", &AindexWrapper::debug_kmer_tf_values,
             "Debug kmer tf values")
        .def("get_tf_values_13mer", &AindexWrapper::get_tf_values_13mer,
             "Get term frequency values for 13-mers")
        .def("get_13mer_tf_array", &AindexWrapper::get_13mer_tf_array,
             "Get direct access to 13-mer tf array")
        .def("get_tf_by_index_13mer", &AindexWrapper::get_tf_by_index_13mer,
             "Get tf value by direct array index for 13-mers")
        .def("get_index_info", &AindexWrapper::get_index_info,
             "Get statistics about loaded index")
        .def("get_positions_13mer", &AindexWrapper::get_positions_13mer,
             "Get positions for 13-mers using the position index");
}
