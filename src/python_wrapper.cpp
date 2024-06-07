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


typedef std::atomic<uint8_t> ATOMIC_BOOL;
emphf::stl_string_adaptor str_adapter;

class UsedReads {
public:
    UsedReads(size_t n_reads) {
        used_reads = new ATOMIC_BOOL[n_reads];
        for (size_t i=0; i < n_reads; ++i) {
            used_reads[i] = 0;
        }
    }

    ~UsedReads() {
        delete[] used_reads;
    }

    ATOMIC_BOOL get(size_t pos) {
        return used_reads[pos].load();
    }  

    void set(size_t pos) {
        ++used_reads[pos];
    }

    bool used_or_use(size_t pos) {
        auto status = get(pos); // 1 is used // atomic
        bool result = true;
        if (status == 0) {
            set(pos);
            result = false;
        }
        return result;
    }

    bool used(size_t pos) {
        auto status = get(pos); // 1 is used // atomic
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
    size_t rid;
    size_t start;
    std::string read;
    size_t pos;
    int ori;
    bool rev;
};

class AindexWrapper {

    size_t *positions = nullptr;
    size_t *indices = nullptr;
    uint32_t *read_pos_cache = nullptr;
    size_t n = 0;
    uint32_t max_tf = 0;
    size_t indices_length = 0;

public:

    bool aindex_loaded = false;
    PHASH_MAP *hash_map;
    size_t n_reads = 0;
    size_t n_kmers = 0;
    std::unordered_map<size_t, uint32_t> start2rid;
    std::vector<size_t> start_positions;

    size_t reads_length = 0;
    char *reads = nullptr;

    std::unordered_map<size_t, size_t> start2end;

    AindexWrapper() {

    }

    ~AindexWrapper() {
        // emphf::logger() << "NOTE: Calling aindex deconstructor..." << std::endl;
        if (positions != nullptr) munmap(positions, n*sizeof(size_t));
        if (indices != nullptr) munmap(indices, indices_length);
        if (reads != nullptr) munmap(reads, reads_length);
        if (read_pos_cache != nullptr) delete[] read_pos_cache;

        delete hash_map;

        reads = nullptr;
        indices = nullptr;
        positions = nullptr;
        read_pos_cache = nullptr;
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

    void _load_reads(const std::string& reads_file) {
        // Memory map reads
        emphf::logger() << "Memory mapping reads file..." << std::endl;
        std::ifstream fout(reads_file, std::ios::in | std::ios::binary);
        fout.seekg(0, std::ios::end);
        size_t length = fout.tellg();
        fout.close();

        int fd = open(reads_file.c_str(), O_RDONLY);
        if (fd == -1) {
            std::cerr << "Failed to open file" << std::endl;
            exit(1);
        }

        char* reads = static_cast<char*>(mmap(nullptr, length, PROT_READ, MAP_PRIVATE, fd, 0));
        close(fd);
        if (reads == MAP_FAILED) {
            std::cerr << "Failed to memory map file" << std::endl;
            exit(2);
        }

        n_reads = 0;
        reads_length = length;

        // Count newlines in parallel
        #pragma omp parallel for reduction(+:n_reads)
        for (size_t i = 0; i < length; ++i) {
            if (reads[i] == '\n') n_reads += 1;
        }
        emphf::logger() << "\tloaded reads: " << n_reads << std::endl;

        emphf::logger() << "\tbuilding start pos index over reads: " << std::endl;

        size_t rid = 0;
        size_t start = 0;
        size_t total = 0;

        uint16_t local_start = 0;
        read_pos_cache = new uint32_t[reads_length];

        // Vector to store intermediate results from threads
        std::vector<std::pair<size_t, size_t>> start_end_pairs;
        std::vector<std::pair<size_t, size_t>> start_rid_pairs;
        std::vector<size_t> local_starts(reads_length, 0);

        #pragma omp parallel
        {
            std::vector<std::pair<size_t, size_t>> local_start_end_pairs;
            std::vector<std::pair<size_t, size_t>> local_start_rid_pairs;

            #pragma omp for
            for (size_t i = 0; i < length; ++i) {
                if (i % 100000000 == 0) {
                    double progress = static_cast<double>(i + 1) / length;
                    printProgressBar(progress);                
                }
                read_pos_cache[i] = local_start;
                if (reads[i] == '\n') {
                    local_start_end_pairs.emplace_back(start, i);
                    local_start_rid_pairs.emplace_back(start, rid);
                    start = i + 1;
                    rid += 1;
                    local_start = 0;
                } else {
                    ++local_start;
                }
            }

            #pragma omp critical
            {
                start_end_pairs.insert(start_end_pairs.end(), local_start_end_pairs.begin(), local_start_end_pairs.end());
                start_rid_pairs.insert(start_rid_pairs.end(), local_start_rid_pairs.begin(), local_start_rid_pairs.end());
            }
        }

        printProgressBar(1.0);

        if (start < length) {
            start_end_pairs.emplace_back(start, total);
            start_rid_pairs.emplace_back(start, rid);
        }

        // Populate maps
        for (const auto& p : start_end_pairs) {
            start2end.emplace(p.first, p.second);
        }
        for (const auto& p : start_rid_pairs) {
            start2rid[p.first] = p.second;
            start_positions.push_back(p.first);
        }

        emphf::logger() << "\tDone" << std::endl;
    }

    void load_reads(std::string reads_file) {
        // Memory map reads
        emphf::logger() << "Memory mapping reads file..." << std::endl;
        std::ifstream fout(reads_file, std::ios::in | std::ios::binary);
        fout.seekg(0, std::ios::end);
        size_t length = fout.tellg();
        fout.close();

        FILE* in = std::fopen(reads_file.c_str(), "rb");
        reads = (char*)mmap(NULL, length, PROT_READ|PROT_WRITE, MAP_PRIVATE, fileno(in), 0);
        if (reads == nullptr) {
            std::cerr << "Failed position loading" << std::endl;
            exit(10);
        }
        fclose(in);

        n_reads = 0;
        reads_length = length;

        for (size_t i=0; i < length; ++i) {
            if (reads[i] == '\n') n_reads += 1;
        }
        emphf::logger() << "\tloaded reads: " << n_reads << std::endl;

        emphf::logger() << "\tbuilding start pos index over reads: " << std::endl;

        size_t rid = 0;
        size_t start = 0;
        size_t total = 0;

        uint16_t local_start = 0;
        read_pos_cache = new uint32_t[reads_length];
        for (size_t i=0; i < length; ++i) {
            if (i % 100000000 == 0 ) {
                double progress = static_cast<double>(i + 1) / length;
                printProgressBar(progress);                
            }
            read_pos_cache[i] = local_start;
            total++;
            if (reads[i] == '\n') {
                start2end.emplace(start, i);
                start2rid[start] = rid;
                start_positions.push_back(start);
                start = i+1;
                rid += 1;
                local_start = 0;
            } else {
                ++local_start;
            }
        }
        printProgressBar(1.0);

        if (start < length) {
            start2end.emplace(start, total);
            start2rid[start] = rid;
            start_positions.push_back(start);
        }

        // for (size_t i=0; i< 500; ++i) {
        //     std::cout << read_pos_cache[i] << " ";
        // }
        std::cout << std::endl;
        
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
        size_t length = fin.tellg();
        fin.seekg(0, std::ios::beg);

        reads = new char[length];
        fin.read(reads, length);
        fin.close();

        if (!reads) {
            std::cerr << "Failed to allocate memory for reads" << std::endl;
            exit(10);
        }

        reads_length = length;
        n_reads = 0;

        for (size_t i = 0; i < length; ++i) {
            if (reads[i] == '\n') n_reads += 1;
        }
        emphf::logger() << "\tloaded reads: " << n_reads << std::endl;

        emphf::logger() << "\tbuilding start pos index over reads: " << std::endl;

        size_t rid = 0;
        size_t start = 0;
        size_t total = 0;

        uint16_t local_start = 0;
        read_pos_cache = new uint32_t[reads_length];
        for (size_t i = 0; i < length; ++i) {
            if (i % 100000000 == 0) {
                double progress = static_cast<double>(i + 1) / length;
                printProgressBar(progress);
            }
            read_pos_cache[i] = local_start;
            total++;
            if (reads[i] == '\n') {
                start2end.emplace(start, i);
                start2rid[start] = rid;
                start_positions.push_back(start);
                start = i + 1;
                rid += 1;
                local_start = 0;
            } else {
                ++local_start;
            }
        }
        printProgressBar(1.0);

        if (start < length) {
            start2end.emplace(start, total);
            start2rid[start] = rid;
            start_positions.push_back(start);
        }

        std::cout << std::endl;
        emphf::logger() << "\tDone" << std::endl;
    }

    void load_aindex(std::string aindex_prefix, uint32_t _max_tf) {
        // Load aindex.

        // std::cout << "Inside load index: " << _max_tf << std::endl;

        n = hash_map->n;
        max_tf = _max_tf;

        std::string pos_file = aindex_prefix + ".pos.bin";
        std::string index_file = aindex_prefix + ".index.bin";
        std::string indices_file = aindex_prefix + ".indices.bin";

        // std::cout << "END" << std::endl;

//        emphf::logger() << "Reading aindex.pos.bin array..." << std::endl;
//
//        size_t f = 0;
//        size_t rid = 0;
//        std::ifstream fout2(pos_file, std::ios::in | std::ios::binary);
//        while(fout2.read(reinterpret_cast<char *>(&f), sizeof(f))) {
//            start_positions.push_back(f);
//            start2rid[f] = rid;
//            rid += 1;
//
//        }
//
//        start_positions.pop_back();
//        start_positions.pop_back();
//        fout2.close();
//        emphf::logger() << "\tDone" << std::endl;

        emphf::logger() << "Reading aindex.indices.bin array..." << std::endl;

        size_t pos = 0;
        std::ifstream fout5(indices_file, std::ios::in | std::ios::binary);
        fout5.seekg(0, std::ios::end);
        size_t length = fout5.tellg();
        fout5.close();

        FILE* in1 = std::fopen(indices_file.c_str(), "rb");
        indices = (size_t*)mmap(NULL, length, PROT_READ|PROT_WRITE, MAP_PRIVATE, fileno(in1), 0);
        if (indices == nullptr) {
            std::cerr << "Failed position loading" << std::endl;
            exit(10);
        }
        fclose(in1);
        indices_length = length;
        emphf::logger() << "\tindices length: " << indices_length << std::endl;
        emphf::logger() << "\tDone" << std::endl;

        emphf::logger() << "Reading aindex.index.bin array..." << std::endl;

        pos = 0;
        std::ifstream fout6(index_file, std::ios::in | std::ios::binary);
        fout6.seekg(0, std::ios::end);
        length = fout6.tellg();
        fout6.close();

        emphf::logger() << "\tpositions length: " << length << std::endl;
        FILE* in = std::fopen(index_file.c_str(), "rb");
        positions = (size_t*)mmap(NULL, length, PROT_READ|PROT_WRITE, MAP_PRIVATE, fileno(in), 0);
        if (positions == nullptr) {
            std::cerr << "Failed position loading" << std::endl;
            exit(10);
        }
        fclose(in);
        this->aindex_loaded = true;
        emphf::logger() << "\tDone" << std::endl;

    }

    std::string get_read_by_rid(uint32_t rid) {
        size_t start = start_positions[rid];
        return get_read_by_start(start);
    }

    size_t get_start_by_pos(size_t pos) {
        return pos - read_pos_cache[pos];
    }

    size_t get_end_by_start(size_t start) {
        return start2end.at(start);
    }

    std::string get_read_by_start(size_t start) const {
        size_t end = start2end.at(start);
        char read[end-start+1];
        std::memcpy(read, &reads[start], end-start);
        read[end-start] = '\0';
        return std::string(read);
    }

    std::string get_read_by_rid(uint32_t rid, int ori) {

        size_t start = start_positions[rid];
        size_t end = start;
        size_t spring_pos = 0;

        while (true) {
            if (reads[end] == '\n') {
                char rkmer[end-spring_pos];
                std::memcpy(rkmer, &reads[spring_pos+1], end-spring_pos-1);
                rkmer[end-spring_pos-1] = '\0';
                return std::string(rkmer);
            } else if (reads[end] == '~') {
                if (ori == 0) {
                    char lkmer[end-start+1];
                    std::memcpy(lkmer, &reads[start], end-start);
                    lkmer[end-start] = '\0';
                    return std::string(lkmer);
                } else {
                    spring_pos = end;
                }
            }
            end += 1;
        }
    }

    size_t get_n() {
        return hash_map->n;
    }

    size_t get_rid(size_t pos) {
        // Get rid by position.
        while (true){
            if (reads[pos] == '\n') {
                return pos+1;
            }
            if (pos == 0) {
                return pos;
            }
            pos -= 1;
        }
    }

    size_t get_start(size_t pos) {
        // Get rid by position.
        uint16_t shift = read_pos_cache[pos];
        return pos - shift;
    }

    size_t get_kid_by_kmer(std::string _kmer) {
        uint64_t kmer = get_dna23_bitset(_kmer);
        return hash_map->get_pfid_by_umer_safe(kmer);
    }

    
    void get_positions(size_t* r, const std::string_view& kmer) {
        // Get read positions and save them to given r
        auto h1 = hash_map->get_pfid(kmer);
        size_t j = 0;
        for (size_t i=indices[h1]; i < indices[h1+1] && h1+1 < indices_length; ++i) {
            if (j == max_tf - 1) {
                break;
            }
            r[j] = positions[i];
            j += 1;
        }
        r[j] = 0;
    }

    size_t get_positions(size_t* r, size_t hash) {
        // Get read positions and save them to given r
        size_t j = 0;
        for (size_t i=indices[hash]; i < indices[hash+1] && hash+1 < indices_length; ++i) {
            if (j == max_tf - 1) {
                break;
            }
            if (positions[i] == 0) { // probably errorprone
                break;
            }
            r[j] = positions[i]-1;
            j += 1;
        }
        r[j] = 0;
        return j;
    }

    std::vector<size_t> get_positions(const std::string& kmer) {
        // Get read positions and save them to given r
        std::vector<size_t> r;
        auto h1 = hash_map->get_pfid(kmer);
        for (size_t i=indices[h1]; i < indices[h1+1] && h1+1 < indices_length; ++i) {
            r.push_back(positions[i]);
        }
        return r;
    }

    std::vector<size_t> get_positions_vector(const std::string& kmer) {
        // Get read positions and save them to given r
        std::vector<size_t> r;
        auto h1 = hash_map->get_pfid(kmer);
        for (size_t i=indices[h1]; i < indices[h1+1] && h1+1 < indices_length; ++i) {
            size_t pos = positions[i];
            if (pos > 0) {
                r.push_back(pos-1);
            }
        }
        return r;
    }

    void get_positions_vector_as_param(const std::string& kmer, std::vector<size_t>& r, size_t max_hits) {
        // Get read positions and save them to given r
        auto h1 = hash_map->get_pfid(kmer);
        size_t total = 0;
        for (size_t i=indices[h1]; i < indices[h1+1] && h1+1 < indices_length; ++i) {
            size_t pos = positions[i];
            if (pos > 0) {
                r.push_back(pos-1);
                total++;
                if (total > max_hits) {
                    break;
                }
            }
        }
    }

    std::unordered_map<size_t, std::vector<size_t> > get_rid2poses(const std::string& kmer) {
        // ''' Wrapper that handle case when two kmer hits in one read.
        // Return rid->poses_in_read dictionary for given kmer. 
        // In this case rid is the start position in reads file.
        // '''
        std::vector<size_t> poses = get_positions_vector(kmer);
        std::unordered_map<size_t, std::vector<size_t> > result;

        for (size_t pos: poses) {
            size_t start = get_start(pos);
            if (!result.count(start)) {
                std::vector<size_t> hits;
                result.emplace(start, hits);
            }
            result[start].push_back(pos - start);
        }
        return result;
    }

    std::vector<Hit> get_reads_by_kmer(const std::string& kmer, UsedReads &used_reads, const size_t min_shift, const size_t max_shift) const {
        uint k = 23;
        std::vector<Hit> result;
        auto h1 = hash_map->get_pfid(kmer);
        for (size_t i=indices[h1]; i < indices[h1+1] && h1+1 < indices_length; ++i) {
            size_t reads_pos = positions[i] - 1;
            uint16_t shift = read_pos_cache[reads_pos]; // zero reserved for empty            
            if (shift > max_shift || shift < min_shift) {
                continue;
            }
            size_t start = reads_pos - shift;
            size_t rid = start2rid.at(start);
            if (used_reads.used(rid)) {
                continue;
            }
            std::string read = get_read_by_start(start);
            std::string check_kmer = read.substr(shift, k);
            bool rev = false;
            if (kmer != check_kmer) {  
                std::string rread = read;
                get_revcomp(read, rread);
                read = rread;
                shift = read.size() - shift - k;
                rev = true;
            }
            Hit hit;
            hit.rid = rid;
            hit.start = start;
            hit.pos = shift;
            hit.rev = rev;
            hit.read = read;
            result.emplace_back(hit);
        }
        return result;
    }

    size_t get(char* ckmer) {
        // Return tf for given char * kmer
        std::string kmer = std::string(ckmer);
        return get(kmer);
    }

    size_t get(uint64_t ukmer) {
        return hash_map->tf_values[ukmer];
    }

    size_t get(std::string kmer) {
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

    size_t get_hash_value(const std::string kmer) {
        // Return hash value for given kmer
        return hash_map->get_pfid(kmer);
    }

     size_t get_hash_value(std::string_view kmer) {
        // Return hash value for given kmer
        return hash_map->get_pfid(kmer);
    }

    size_t get(std::string_view kmer) const {
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

    size_t get_strand(const std::string& kmer) {
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
                return 1;
            }
        } else {
            return 2;
        }
        return 0;
    }

    void get_kmer_by_kid(size_t r, char* kmer) {
            // if (r >= hash_map->n) {
            //     return;
            // }
            uint64_t ukmer = hash_map->checker[r];
            get_bitset_dna23_c(ukmer, kmer, 23);            
    }

    size_t get_kmer(size_t p, char* kmer, char* rkmer) {
        // Get tf, kmer and rev_kmer stored in given arrays.
        uint64_t ukmer = hash_map->checker[p];
        uint64_t urev_kmer = reverseDNA(ukmer);
        get_bitset_dna23_c(ukmer, kmer, 23);
        get_bitset_dna23_c(urev_kmer, rkmer, 23);
        return hash_map->tf_values[p];
    }


    size_t get_hash_size() {
        return hash_map->n;
    }

    void increase(char* ckmer) {
        std::string kmer = std::string(ckmer);
        hash_map->increase(kmer);
    }

    void decrease(char* ckmer) {
        std::string kmer = std::string(ckmer);
        hash_map->decrease(kmer);
    }

    void set_positions(size_t* r, const std::string& kmer) {
        // Set read positions
        auto h1 = hash_map->get_pfid(kmer);
        size_t j = 0;
        for (size_t i=indices[h1]; i < indices[h1+1]; ++i) {
            positions[i] = r[j];
            j += 1;
        }
    }

    void check_aindex() {

        for (size_t h1=0; h1<hash_map->n; ++h1) {
            size_t tf = hash_map->tf_values[h1];
            size_t xtf = 0;

            if (h1 && h1 % 1000000 == 0) {
                std::cout << "Completed: " << h1 << "/" << hash_map->n << std::endl;
            }

            for (size_t i=indices[h1]; i < indices[h1+1]; ++i) {
                if (positions[i] == 0) {
                    break;
                }

                xtf += 1;

                size_t pos = positions[i]-1;

                char ckmer[Settings::K];

                std::memcpy(ckmer, &reads[pos], Settings::K);
                ckmer[Settings::K] = '\0';
                std::string data_kmer = std::string(ckmer);

                uint64_t h1_kmer = hash_map->checker[h1];
                std::string kmer = get_bitset_dna23(h1_kmer);
                if (data_kmer != kmer) {
                    size_t rh1 = reverseDNA(h1_kmer);
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

        for (size_t h1=0; h1<hash_map->n; ++h1) {

            if (h1 && h1 % 1000000 == 0) {
                std::cout << "Completed: " << h1 << "/" << hash_map->n << std::endl;
            }

            uint64_t h1_kmer = hash_map->checker[h1];
            std::string kmer = get_bitset_dna23(h1_kmer);
            hits.clear();
            get_reads_se_by_kmer(kmer, h1, used_reads, hits);

            size_t max_pos = 0;

            for (auto hit: hits) {
                std::max(max_pos, hit.pos);
                std::string subkmer = hit.read.substr(hit.pos, Settings::K);
                assert(subkmer == kmer);
                std::cout << kmer << " " << subkmer << " " << h1 << " " << hash_map->tf_values[h1] << std::endl;
            }
        }
    }


    void get_reads_se_by_kmer(std::string const kmer, size_t h1, bool* used_reads, std::vector<Hit> &hits) {


        for (size_t i=indices[h1]; i < indices[h1+1]; ++i) {

            if (positions[i] == 0) {
                break;
            }

            size_t position = positions[i] - 1;
            size_t start = get_start(position);

            size_t end = start;
            size_t spring_pos = 0;

            size_t pos = position - start;
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

            size_t real_rid = start2rid[start];

            Hit hit;
            hit.rid = real_rid;
            hit.start = start;
            hit.pos = pos;
            hit.rev = 0;
            spring_pos = spring_pos - start;

            if (pos < spring_pos) {
                hit.read = left_read;
                hit.ori = 0;

                if (hit.read.substr(hit.pos, Settings::K) != kmer) {
                    std::string rleft_read = hit.read;
                    get_revcomp(hit.read, rleft_read);
                    hit.pos = hit.read.length() - pos - Settings::K;
                    if (rleft_read.substr(hit.pos, Settings::K) != kmer) {
                        std::cout << rleft_read << std::endl;
                        std::cout << left_read << std::endl;
                        std::cout << right_read << std::endl;
                        std::cout << kmer << " " << pos <<  std::endl;
                        continue;
                    }
                    hit.read = rleft_read;
                    hit.rev = 1;
//                    hit.ori = 1;
                }
            } else {

                if (hit.pos == spring_pos) {
                    hit.pos = hit.pos - spring_pos;
                    std::cout <<  left_read << std::endl;
                    std::cout <<  right_read << std::endl;
                    std::cout << kmer << std::endl;

                } else {
                    hit.pos = hit.pos - spring_pos - 1;
                }

                hit.read = right_read;
                hit.ori = 1;

                if (hit.read.substr(hit.pos, Settings::K) != kmer) {
                    std::string rright_read = hit.read;
                    get_revcomp(hit.read, rright_read);
                    hit.pos = hit.read.length() - hit.pos - Settings::K;
                    if (rright_read.substr(hit.pos, Settings::K) != kmer) {
                        continue;
                    }
                    hit.read = rright_read;
                    hit.rev = 1;
//                    hit.ori = 0;
                }

            }

//            if (reversed_reads[real_rid]) {
//                if (hit.ori) {
//                    hit.ori = 0;
//                } else {
//                    hit.ori = 1;
//                }
//            }

            if (used_reads[2*hit.rid+hit.ori]) {
                continue;
            }
            hits.push_back(hit);

        }
    }


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
                const size_t max_tf,
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

    void AindexWrapper_load_index(AindexWrapper* foo, char* index_prefix, uint32_t max_tf){ foo->load_aindex(index_prefix, max_tf); }
    
    void AindexWrapper_increase(AindexWrapper* foo, char* kmer){ foo->increase(kmer); }
    void AindexWrapper_decrease(AindexWrapper* foo, char* kmer){ foo->decrease(kmer); }

    size_t AindexWrapper_get_kid_by_kmer(AindexWrapper* foo, char* kmer){ return foo->get_kid_by_kmer(kmer); }

    void AindexWrapper_get_kmer_by_kid(AindexWrapper* foo, size_t kid, char* kmer){ foo->get_kmer_by_kid(kid, kmer); }

    size_t AindexWrapper_get(AindexWrapper* foo, char* kmer){ return foo->get(kmer); }

    size_t AindexWrapper_get_n(AindexWrapper* foo){ return foo->get_n(); }

    size_t AindexWrapper_get_rid(AindexWrapper* foo, size_t pos){ return foo->get_rid(pos); }

//    char* AindexWrapper_get_read(AindexWrapper* foo, size_t start, int ori){ return foo->get_read(pos, ori); }

    void AindexWrapper_get_positions(AindexWrapper* foo, size_t* r, char* kmer){ foo->get_positions(r, kmer); }


    void AindexWrapper_set_positions(AindexWrapper* foo, size_t* r, char* kmer){ foo->set_positions(r, kmer); }

    size_t AindexWrapper_get_kmer(AindexWrapper* foo, size_t p, char* kmer, char* rkmer){ return foo->get_kmer(p, kmer, rkmer); }

    size_t AindexWrapper_get_strand(AindexWrapper* foo, char* kmer){ return foo->get_strand(kmer); }

    size_t AindexWrapper_get_hash_size(AindexWrapper* foo){ return foo->get_hash_size(); }
}

#endif