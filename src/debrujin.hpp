//
// Created by Aleksey Komissarov on 30/08/15.
//
#ifndef STIRKA_DEBRUJIN_H
#define STIRKA_DEBRUJIN_H

#include <stdint.h>
#include "hash.hpp"
#include "read.hpp"
#include <sstream>

namespace DEBRUJIN {

    struct CONT {
        uint32_t A = 0;
        uint32_t C = 0;
        uint32_t G = 0;
        uint32_t T = 0;
        uint32_t n = 0;
        uint32_t sum = 0;
        char best_hit;
        uint32_t best_hit_tf = 0;
        uint64_t best_ukmer = 0;

        void print(int i) {
            std::cout << "POS: " << i << " N: " << n << " A: " << A << " C: " << C << " G: " << G << " T: " << T << std::endl;
        }

        std::string get_cont() {
            std::stringstream ss;
            ss << A << " " << C << " " << G << " " << T;
            return ss.str();
        }
    };


    int get_freq(uint64_t kmer, PHASH_MAP &kmers);

    void print_next(uint64_t kmer, PHASH_MAP &kmers, CONT &cont, uint32_t cutoff);

    void print_prev(uint64_t kmer, PHASH_MAP &kmers, CONT &cont, uint32_t cutoff);

    void set_fm_for_read(READS::READ &read, PHASH_MAP &kmers);

    void set_fm_for_read(READS::READ &read, PHASH_MAP &kmers, int from, uint64_t to);


}
#endif //STIRKA_DEBRUJIN_H
