//
// Created by Aleksey Komissarov on 15/01/2020.
//

int main(int argc, char** argv) {

    if (argc < 5) {
        std::cerr << "Build position index over reference with annotation." << std::endl;
        std::cerr << "Expected arguments: " << argv[0]
                  << " reference.reads reference.header pf_prefix output_prefix" << std::endl;
        std::terminate();
    }

    std::string reads_file = argv[1];
    std::string header file = argv[2];
    std::string pf_prefix = argv[3];
    std::string output_prefix = argv[4];

    std::string pf_file = pf_prefix + ".pf";
    std::string kmer_bin_file = pf_prefix + ".kmers.bin";

    std::string output_tf_bin;
    std::string output_positions_bin;
    std::string output_header_txt;

    PHASH_MAP hash_map;

    // 0. Load hash

    load_hash_with_empty_tf(hash_map, pf_prefix, pf_prefix);

    // 1. Count tf

    count_kmer_in_reads_in_parallel(hash_map, reads_file);

    // 2. Init: (refid, pos) data structure

    // 3. Fill it with data over reads

    // 4. Save data

    // Usage: kmer => chr1 20, chr2 45 and TF=2

}