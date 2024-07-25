#ifndef AINDEX_FILE_H
#define AINDEX_FILE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/mman.h>
#include <atomic>
#include <mutex>
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
    void addInterval(uint64_t rid, uint64_t start, uint64_t end) {
        intervals.insert({rid, start, end});
    }

    // Finds the interval that contains the position
    const Interval* findInterval(uint64_t pos) const {
        auto it = intervals.lower_bound({0, pos, pos});
        if (it != intervals.begin()) {
            --it;
            if (it->start <= pos && pos <= it->end) {
                return &(*it);
            }
        }
        return nullptr;
    }

private:
    std::set<Interval> intervals;
};

class UsedReads {
public:

    UsedReads(uint64_t n_reads) {
        used_reads = new ATOMIC_BOOL[n_reads];
        for (uint64_t i=0; i < n_reads; ++i) {
            used_reads[i] = 0;
        }
    }

    ~UsedReads() {
        delete[] used_reads;
    }

    ATOMIC_BOOL get(uint64_t rid) {
        return used_reads[rid].load();
    }  

    void set(uint64_t rid) {
        used_reads[rid].store(true);
    }

    bool used_or_use(uint64_t rid) {
        auto status = get(rid); // 1 is used // atomic
        bool result = true;
        if (status == 0) {
            set(rid);
            result = false;
        }
        return result;
    }

    bool used(uint64_t rid) {
        auto status = get(rid); // 1 is used // atomic
        bool result = true;
        if (status == 0) {
            result = false;
        }
        return result;
    }

private:
    ATOMIC_BOOL* used_reads;
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
    
    AindexWrapper() {

    }

    ~AindexWrapper() {
        // emphf::logger() << "NOTE: Calling aindex deconstructor..." << std::endl;
        if (positions != nullptr) munmap(positions, n*sizeof(uint64_t));
        if (indices != nullptr) munmap(indices, indices_length);
        if (reads != nullptr) munmap(reads, reads_size);

        delete hash_map;

        reads = nullptr;
        indices = nullptr;
        positions = nullptr;
    }

    void load(std::string index_prefix, std::string tf_file){

        hash_map = new PHASH_MAP();
        // Load perfect hash into hash_map into memory
        emphf::logger() << "Reading index and hash..." << std::endl;
        std::string hash_filename = index_prefix + ".pf";
        emphf::logger() << "...files: " << index_prefix << std::endl;
        emphf::logger() << "...files: " << tf_file << std::endl;
        emphf::logger() << "...files: " << hash_filename << std::endl;
        load_hash(*hash_map, index_prefix, tf_file, hash_filename);
        n_kmers = hash_map->n;
        emphf::logger() << "\tDone" << std::endl;
    }

    void load_hash_file(std::string hash_filename) {
        emphf::logger() << "Loading only hash..." << std::endl;
        load_only_hash(*hash_map, hash_filename);
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

        if (!reads) {
            std::cerr << "Failed to allocate memory for reads" << std::endl;
            exit(10);
        }

        reads_size = length;

        emphf::logger() << "\tbuilding start pos index over reads: " << std::endl;
        std::string index_file = reads_file.substr(0, reads_file.find_last_of(".")) + ".ridx";
        load_reads_index(index_file);
        emphf::logger() << "\tDone" << std::endl;
    }

    void load_aindex(std::string aindex_prefix, uint32_t _max_tf) {
        // Load aindex.

        n = hash_map->n;
        max_tf = _max_tf;

        std::string pos_file = aindex_prefix + ".pos.bin";
        std::string index_file = aindex_prefix + ".index.bin";
        std::string indices_file = aindex_prefix + ".indices.bin";

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

    uint64_t get_reads_size() {
        return reads_size;
    }

    uint64_t get_n() {
        return hash_map->n;
    }

    uint64_t get_hash_size() {
        return hash_map->n;
    }

    // Various getters for reads

    const char* get_read(uint64_t start, uint64_t end, uint rev) {
        // TODO: make it thread safe
        if (start >= reads_size || end > reads_size || start >= end) {
            return nullptr;  // Invalid range
        }
        static std::string read_str;
        read_str = std::string(reads + start, end - start);
        if (rev > 0) {
            read_str = get_revcomp(read_str);
        }
        return read_str.c_str();
    }

    std::string get_read_by_rid(uint32_t rid) {
        uint64_t start = start_positions[rid];
        uint64_t end = start2end.at(start);
        return std::string(reads + start, end - start);
    }

    const char * get_pointer_to_read_by_rid(uint64_t rid) {
        // TODO: make it thread safe
        if (rid >= start_positions.size()) {
            std::cerr << "Read id " << rid << " not found." << std::endl;
            std::terminate();
        }
        uint64_t start = start_positions[rid];
        uint64_t end = start2end.at(start);
        static std::string read_str;
        read_str = std::string(reads + start, end - start);
        return read_str.c_str();
    }

    uint64_t get_start_by_pos(uint64_t pos) {
        const Interval* interval = pos_intervalTree.findInterval(pos);
        if (interval) {
            return interval->start;
        } else {
            std::cerr << "Position " << pos << " not found in any interval." << std::endl;
            std::terminate();
        }
    }

    uint64_t get_end_by_start(uint64_t start) {
        return start2end.at(start);
    }

    std::string get_read_by_start(uint64_t start) const {
        uint64_t end = start2end.at(start);
        return std::string(reads + start, end - start);;
    }

    uint64_t get_rid(uint64_t pos) {
        const Interval* interval = pos_intervalTree.findInterval(pos);
        if (!interval) {
            std::cerr << "Position " << pos << " not found in any interval." << std::endl;
            std::terminate();
        }
        uint64_t start = interval->start;
        return start2rid.at(start);
    }

    // Varios getters for kmers

    uint64_t get(char* ckmer) {
        // Return tf for given char * kmer
        return get(std::string(ckmer));
    }

    uint64_t get(uint64_t ukmer) {
        if (ukmer >= hash_map->n) {
            return 0;
        }
        return hash_map->tf_values[ukmer];
    }

    uint64_t get(std::string& kmer) {
        // Return tf for given kmer
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

    uint64_t get(std::string_view kmer) const {
        // Return tf for given kmer
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

    uint64_t get_hash_value(std::string_view kmer) {
        // Return hash value for given kmer
        return hash_map->get_pfid(kmer);
    }

    uint64_t get_strand(const std::string& kmer) {
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
                return 2;
            }
        } else {
            return 1;
        }
        return 0;
    }

    void get_kmer_by_kid(uint64_t r, char* kmer) {
            uint64_t ukmer = hash_map->checker[r];
            get_bitset_dna23_c(ukmer, kmer, 23);            
    }

    uint64_t get_kmer(uint64_t kid, char* kmer, char* rkmer) {
        // Get tf, kmer and rev_kmer stored in given arrays.
        // TODO: fix this
        uint64_t ukmer = hash_map->checker[kid];
        uint64_t urev_kmer = reverseDNA(ukmer);
        get_bitset_dna23_c(ukmer, kmer, 23);
        get_bitset_dna23_c(urev_kmer, rkmer, 23);
        return hash_map->tf_values[kid];
    }

    uint64_t get_kid_by_kmer(std::string _kmer) {
        uint64_t kmer = get_dna23_bitset(_kmer);
        return hash_map->get_pfid_by_umer_safe(kmer);
    }
    
    // Getters for positions

    void get_positions(uint64_t* r, const std::string_view& kmer) {
        // Get read positions and save them to given r
        auto h1 = hash_map->get_pfid(kmer);
        uint64_t j = 0;
        for (uint64_t i=indices[h1]; i < indices[h1+1] && h1+1 < indices_length; ++i) {
            if (j == max_tf - 1) {
                break;
            }
            r[j] = positions[i];
            j += 1;
        }
        r[j] = 0;
    }

    std::vector<uint64_t> get_positions(const std::string& kmer) {
        // Get read positions and save them to given r
        std::vector<uint64_t> r;
        auto h1 = hash_map->get_pfid(kmer);
        for (uint64_t i=indices[h1]; i < indices[h1+1] && h1+1 < indices_length; ++i) {
            if (positions[i] == 0) {
                continue;
            }
            r.push_back(positions[i]-1);
        }
        return r;
    }

    // Aindex manipulation

    void increase(char* ckmer) {
        std::string kmer = std::string(ckmer);
        hash_map->increase(kmer);
    }

    void decrease(char* ckmer) {
        std::string kmer = std::string(ckmer);
        hash_map->decrease(kmer);
    }

    void set_positions(uint64_t* r, const std::string& kmer) {
        // Set read positions
        // TODO: check borders
        auto h1 = hash_map->get_pfid(kmer);
        uint64_t j = 0;
        for (uint64_t i=indices[h1]; i < indices[h1+1]; ++i) {
            positions[i] = r[j];
            j += 1;
        }
    }

    // Consistency checks

    void check_get_reads_se_by_kmer(std::string const kmer, uint64_t h1, bool* used_reads, std::vector<Hit> &hits) {

        for (uint64_t i=indices[h1]; i < indices[h1+1]; ++i) {

            if (positions[i] == 0) {
                break;
            }

            uint64_t position = positions[i] - 1;
            uint64_t start = get_start_by_pos(position);

            uint64_t end = start;
            uint64_t spring_pos = 0;

            uint64_t pos = position - start;
            std::string left_read;
            std::string right_read;

            while (true) {
                if (reads[end] == '\n') {
                    if (spring_pos > 0) {
                        char rkmer[end-spring_pos];
                        std::memcpy(rkmer, &reads[spring_pos+1], end-spring_pos-1);
                        rkmer[end-spring_pos-1] = '\0';
                        right_read = std::string(rkmer);
                    }
                    break;
                } else if (reads[end] == '~') {
                    char lkmer[end-start+1];
                    std::memcpy(lkmer, &reads[start], end-start);
                    lkmer[end-start] = '\0';
                    left_read = std::string(lkmer);
                    spring_pos = end;
                }
                end += 1;
            }

            uint64_t real_rid = start2rid[start];

            Hit hit;
            hit.rid = real_rid;
            hit.start = start;
            hit.local_pos = pos;
            hit.rev = 0;
            spring_pos = spring_pos - start;

            if (pos < spring_pos) {
                hit.read = left_read;
                hit.ori = 0;

                if (hit.read.substr(hit.local_pos, Settings::K) != kmer) {
                    std::string rleft_read = hit.read;
                    get_revcomp(hit.read, rleft_read);
                    hit.local_pos = hit.read.length() - pos - Settings::K;
                    if (rleft_read.substr(hit.local_pos, Settings::K) != kmer) {
                        std::cout << rleft_read << std::endl;
                        std::cout << left_read << std::endl;
                        std::cout << right_read << std::endl;
                        std::cout << kmer << " " << pos <<  std::endl;
                        continue;
                    }
                    hit.read = rleft_read;
                    hit.rev = 1;
                }
            } else {

                if (hit.local_pos == spring_pos) {
                    hit.local_pos = hit.local_pos - spring_pos;
                    std::cout <<  left_read << std::endl;
                    std::cout <<  right_read << std::endl;
                    std::cout << kmer << std::endl;

                } else {
                    hit.local_pos = hit.local_pos - spring_pos - 1;
                }

                hit.read = right_read;
                hit.ori = 1;

                if (hit.read.substr(hit.local_pos, Settings::K) != kmer) {
                    std::string rright_read = hit.read;
                    get_revcomp(hit.read, rright_read);
                    hit.local_pos = hit.read.length() - hit.local_pos - Settings::K;
                    if (rright_read.substr(hit.local_pos, Settings::K) != kmer) {
                        continue;
                    }
                    hit.read = rright_read;
                    hit.rev = 1;
                }

            }

            if (used_reads[2*hit.rid+hit.ori]) {
                continue;
            }
            hits.push_back(hit);

        }
    }

    void check_aindex() {

        for (uint64_t h1=0; h1<hash_map->n; ++h1) {
            uint64_t tf = hash_map->tf_values[h1];
            uint64_t xtf = 0;

            if (h1 && h1 % 1000000 == 0) {
                std::cout << "Completed: " << h1 << "/" << hash_map->n << std::endl;
            }

            for (uint64_t i=indices[h1]; i < indices[h1+1]; ++i) {
                if (positions[i] == 0) {
                    break;
                }

                xtf += 1;

                uint64_t pos = positions[i]-1;

                char ckmer[Settings::K];

                std::memcpy(ckmer, &reads[pos], Settings::K);
                ckmer[Settings::K] = '\0';
                std::string data_kmer = std::string(ckmer);

                uint64_t h1_kmer = hash_map->checker[h1];
                std::string kmer = get_bitset_dna23(h1_kmer);
                if (data_kmer != kmer) {
                    uint64_t rh1 = reverseDNA(h1_kmer);
                    std::string rkmer = get_bitset_dna23(rh1);
                    if (data_kmer != rkmer) {
                        std::cout << h1 << " " << i << " " << tf << " " << xtf << " " <<  data_kmer << " " << kmer << " " << rkmer << std::endl;
                    }
                }
            }
            if (tf != xtf) {
                std::cout << tf << " " << xtf << std::endl;

            }
        }
    }

    void check_aindex_reads() {

        bool* used_reads = new bool[1];
        std::vector<Hit> hits;

        for (uint64_t h1=0; h1<hash_map->n; ++h1) {

            if (h1 && h1 % 1000000 == 0) {
                std::cout << "Completed: " << h1 << "/" << hash_map->n << std::endl;
            }

            uint64_t h1_kmer = hash_map->checker[h1];
            std::string kmer = get_bitset_dna23(h1_kmer);
            hits.clear();
            check_get_reads_se_by_kmer(kmer, h1, used_reads, hits);

            uint64_t max_pos = 0;

            for (auto hit: hits) {
                std::max(max_pos, hit.local_pos);
                std::string subkmer = hit.read.substr(hit.local_pos, Settings::K);
                assert(subkmer == kmer);
                std::cout << kmer << " " << subkmer << " " << h1 << " " << hash_map->tf_values[h1] << std::endl;
            }
        }
    }

    // Deconstructors

    void freeme(char* ptr)
    {
        std::cout << "freeing address: " << ptr << std::endl;
        free(ptr);
    }
};

AindexWrapper load_aindex(
                const std::string index_prefix,
                const std::string tf_prefix,
                const std::string input_reads_file,
                const std::string aindex_prefix,
                const uint64_t max_tf,
                bool in_memory = false
                ) {
    AindexWrapper aindex = AindexWrapper();
    std::string tf_file = tf_prefix + ".tf.bin";
    aindex.load(index_prefix, tf_file);
    if (in_memory) {
        aindex.load_reads_in_memory(input_reads_file);
    } else {
        aindex.load_reads(input_reads_file);
    }
    aindex.load_aindex(aindex_prefix, max_tf);
    return aindex;
}

AindexWrapper load_index(
                const std::string index_prefix,
                const std::string tf_prefix
                ) {
    AindexWrapper aindex = AindexWrapper();
    std::string tf_file = tf_prefix + ".tf.bin";
    aindex.load(index_prefix, tf_file);
    return aindex;
}

extern "C" {

    AindexWrapper* AindexWrapper_new(){ return new AindexWrapper(); }

    void AindexWrapper_load(AindexWrapper* foo, char* index_prefix, char* tf_file){ foo->load(index_prefix, tf_file); }

    void AindexWrapper_freeme(AindexWrapper* foo, char* ptr){ foo->freeme(ptr); }

    void AindexWrapper_load_hash_file(AindexWrapper* foo, char* hash_filename, char* tf_file){ foo->load(hash_filename, tf_file); }

    void AindexWrapper_load_reads(AindexWrapper* foo, char* reads_file){ foo->load_reads(reads_file); }

    void AindexWrapper_load_reads_index(AindexWrapper* foo, char* index_file){ foo->load_reads_index(index_file); }

    void AindexWrapper_load_index(AindexWrapper* foo, char* index_prefix, uint32_t max_tf){ foo->load_aindex(index_prefix, max_tf); }
    
    void AindexWrapper_increase(AindexWrapper* foo, char* kmer){ foo->increase(kmer); }

    void AindexWrapper_decrease(AindexWrapper* foo, char* kmer){ foo->decrease(kmer); }

    uint64_t AindexWrapper_get_kid_by_kmer(AindexWrapper* foo, char* kmer){ return foo->get_kid_by_kmer(kmer); }

    void AindexWrapper_get_kmer_by_kid(AindexWrapper* foo, uint64_t kid, char* kmer){ foo->get_kmer_by_kid(kid, kmer); }

    uint64_t AindexWrapper_get(AindexWrapper* foo, char* kmer){ return foo->get(kmer); }

    uint64_t AindexWrapper_get_n(AindexWrapper* foo){ return foo->get_n(); }

    uint64_t AindexWrapper_get_rid(AindexWrapper* foo, uint64_t pos){ return foo->get_rid(pos); }

    uint64_t AindexWrapper_get_start(AindexWrapper* foo, uint64_t pos){ return foo->get_start_by_pos(pos); }

    const char*  AindexWrapper_get_read(AindexWrapper* foo, uint64_t start, uint64_t end, uint rev){ return foo->get_read(start, end, rev); }

    const char*  AindexWrapper_get_read_by_rid(AindexWrapper* foo, uint64_t rid){ return foo->get_pointer_to_read_by_rid(rid); }

    void AindexWrapper_get_positions(AindexWrapper* foo, uint64_t* r, char* kmer){ foo->get_positions(r, kmer); }

    void AindexWrapper_set_positions(AindexWrapper* foo, uint64_t* r, char* kmer){ foo->set_positions(r, kmer); }

    uint64_t AindexWrapper_get_kmer(AindexWrapper* foo, uint64_t p, char* kmer, char* rkmer){ return foo->get_kmer(p, kmer, rkmer); }

    uint64_t AindexWrapper_get_strand(AindexWrapper* foo, char* kmer){ return foo->get_strand(kmer); }

    uint64_t AindexWrapper_get_hash_size(AindexWrapper* foo){ return foo->get_hash_size(); }

    uint64_t AindexWrapper_get_reads_size(AindexWrapper* foo){ return foo->get_reads_size(); }

    void AindexWrapper_load_reads_in_memory(AindexWrapper* foo, char* reads_file){ foo->load_reads_in_memory(reads_file); }

    void AindexWrapper_load_aindex(AindexWrapper* foo, char* aindex_prefix, uint32_t max_tf){ foo->load_aindex(aindex_prefix, max_tf); }

    uint64_t AindexWrapper_get_start_by_pos(AindexWrapper* foo, uint64_t pos) { return foo->get_start_by_pos(pos); }

    uint64_t AindexWrapper_get_end_by_start(AindexWrapper* foo, uint64_t start) { return foo->get_end_by_start(start); }

    const char* AindexWrapper_get_pointer_to_read_by_rid(AindexWrapper* foo, uint64_t rid) { return foo->get_pointer_to_read_by_rid(rid); }

    uint64_t AindexWrapper_get_hash_value(AindexWrapper* foo, char* kmer) { return foo->get_hash_value(kmer); }

    void AindexWrapper_check_aindex(AindexWrapper* foo) { foo->check_aindex(); }

    void AindexWrapper_check_aindex_reads(AindexWrapper* foo) { foo->check_aindex_reads(); }
}

#endif