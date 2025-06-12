#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <thread>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits.h>
#include <math.h>
#include <mutex>
#include <string_view>
#include <filesystem>
#include "emphf/common.hpp"
#include "read.hpp"

int main(int argc, char** argv) {

    if (argc < 5) {
        std::cerr << "Convert fasta or fastq reads to simple reads." << std::endl;
        std::cerr << "Expected arguments: " << argv[0]
        << " <fastq_file1|fasta_file1|reads_file> <fastq_file2|-> <fastq|fasta|se|reads> <output_prefix>" << std::endl;
        std::terminate();
    }

    std::string file_name1 = argv[1];
    std::string file_name2 = argv[2];
    std::string read_type = argv[3];
    std::string output_prefix = argv[4];
    
    // Check and create output directory if needed
    std::filesystem::path output_path(output_prefix);
    std::filesystem::path output_dir = output_path.parent_path();
    
    if (!output_dir.empty() && !std::filesystem::exists(output_dir)) {
        emphf::logger() << "Output directory does not exist. Creating: " << output_dir << std::endl;
        try {
            std::filesystem::create_directories(output_dir);
            emphf::logger() << "Successfully created directory: " << output_dir << std::endl;
        } catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error creating directory " << output_dir << ": " << e.what() << std::endl;
            std::terminate();
        }
    }
    
    std::string index_file = output_prefix + ".ridx";
    std::string header_file = output_prefix + ".header";
    std::string output_file = output_prefix + ".reads";

    emphf::logger() << "Starting..." << std::endl;

    emphf::logger() << "Converting reads..." << std::endl;
    uint64_t n_reads = 0;

    std::string line1;
    std::string line2;

    if (read_type == "fastq") {
        std::ofstream fout(output_file, std::ios::out);
        std::ofstream fout_index(index_file, std::ios::out);
        std::ifstream fin1(file_name1, std::ios::in);
        std::ifstream fin2(file_name2, std::ios::in);

        uint64_t start_pos = 0;
        while (std::getline(fin1, line1)) {
            std::getline(fin1, line1);
            std::getline(fin2, line2);
            std::getline(fin2, line2);

            uint64_t end_pos = start_pos + line1.size() + line2.size() + 1; // Adding 1 for the '~' separator

            std::string revcomp_line2 = get_revcomp(line2);

            fout << line1;
            fout << "~";
            fout << revcomp_line2;
            fout << "\n";

            fout_index << n_reads << "\t" << start_pos << "\t" << end_pos << "\n";

            start_pos = end_pos + 1; // Adding 1 for the newline character

            std::getline(fin1, line1);
            std::getline(fin1, line1);
            std::getline(fin2, line2);
            std::getline(fin2, line2);
            n_reads += 1;

            if (n_reads % 1000000 == 0) {
                emphf::logger() << "Completed: " << n_reads << std::endl;
            }
        }

        fin1.close();
        fin2.close();
        fout.close();
        fout_index.close();

    } else if (read_type == "se") {
        std::ofstream fout(output_file, std::ios::out);
        std::ofstream fout_index(index_file, std::ios::out);
        std::ifstream fin1(file_name1, std::ios::in);

        uint64_t start_pos = 0;
        while (std::getline(fin1, line1)) {
            std::getline(fin1, line1);
            
            uint64_t end_pos = start_pos + line1.size();

            fout << line1;
            fout << "\n";

            fout_index << n_reads << "\t" << start_pos << "\t" << end_pos << "\n";

            start_pos = end_pos + 1; // Adding 1 for the newline character

            std::getline(fin1, line1);
            std::getline(fin1, line1);
            n_reads += 1;

            if (n_reads % 1000000 == 0) {
                emphf::logger() << "Completed: " << n_reads << std::endl;
            }
        }

        fin1.close();
        fout.close();
        fout_index.close();

    } else if (read_type == "reads") {
        std::ifstream fin1(file_name1, std::ios::in);
        std::ofstream fout_index(index_file, std::ios::out);
        uint64_t start_pos = 0;
        while (std::getline(fin1, line1)) {
            
            uint64_t end_pos = start_pos + line1.size();
            fout_index << n_reads << "\t" << start_pos << "\t" << end_pos << "\n";

            start_pos = end_pos + 1; // Adding 1 for the newline character

            n_reads += 1;

            if (n_reads % 1000000 == 0) {
                emphf::logger() << "Completed: " << n_reads << std::endl;
            }
        }

        fin1.close();
        fout_index.close();
        
    } else if (read_type == "fasta") {
        std::ofstream fout(output_file, std::ios::out);
        std::ofstream fout_index(index_file, std::ios::out);
        std::ifstream fin1(file_name1, std::ios::in);

        std::ofstream fout_header(header_file, std::ios::out);

        std::string current_sequence;
        std::string header;
        uint64_t start_pos = 0;
        while (std::getline(fin1, line1)) {
            if (line1[0] == '>') {
                if (!current_sequence.empty()) {
                    uint64_t end_pos = start_pos + current_sequence.size();

                    fout << current_sequence << "\n";
                    fout_index << n_reads << "\t" << start_pos << "\t" << end_pos << "\n";
                    fout_header << header << "\t" << start_pos << "\t" <<  current_sequence.size() << "\n";

                    start_pos = end_pos + 1; // Adding 1 for the newline character
                    n_reads += 1;
                    current_sequence.clear();

                    if (n_reads % 1000000 == 0) {
                        emphf::logger() << "Completed: " << n_reads << std::endl;
                    }
                }
                header = line1.substr(1);
                continue;
            }
            current_sequence += line1;
        }

        if (!current_sequence.empty()) {
            uint64_t end_pos = start_pos + current_sequence.size();

            fout << current_sequence << "\n";
            fout_index << n_reads << "\t" << start_pos << "\t" << end_pos << "\n";
            fout_header << header << "\t" << start_pos << "\t" <<  current_sequence.size() << "\n";
            
            n_reads += 1;
        }

        fin1.close();
        fout_header.close();
        fout.close();
        fout_index.close();
    } else {
        emphf::logger() << "Unknown format." << std::endl;
        exit(2);
    }

    

    return 0;
}