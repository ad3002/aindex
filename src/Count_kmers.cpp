#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <algorithm>
#include <memory>
#include <cstring>
#include <immintrin.h>
#include <iomanip>

// Версия с локальными хеш-таблицами для каждого потока
using kmer_t = uint64_t;

// Оптимизированная хеш-функция для k-меров
struct KmerHash {
    size_t operator()(kmer_t k) const {
        // MurmurHash-inspired mixing
        k ^= k >> 33;
        k *= 0xff51afd7ed558ccdULL;
        k ^= k >> 33;
        k *= 0xc4ceb9fe1a85ec53ULL;
        k ^= k >> 33;
        return k;
    }
};

// Используем robin_hood map или просто резервируем больше места
template<typename K, typename V>
using FastMap = std::unordered_map<K, V, KmerHash>;

class OptimizedKmerCounter {
private:
    size_t k;
    size_t num_threads;
    size_t min_count;
    bool use_canonical;
    static constexpr size_t BATCH_SIZE = 1000;  // Обрабатываем пачками
    static constexpr size_t BUFFER_SIZE = 16 * 1024 * 1024; // 16MB буфер для чтения
    
    // Каждый поток имеет свою локальную таблицу
    struct ThreadData {
        FastMap<kmer_t, uint32_t> local_kmers;
        std::vector<std::string> batch;
        size_t sequences_processed = 0;
        size_t kmers_processed = 0;
        char* read_buffer = nullptr;
        
        ThreadData() {
            local_kmers.reserve(1 << 20); // Резервируем место для 1M k-меров
            batch.reserve(BATCH_SIZE);
            read_buffer = new char[BUFFER_SIZE];
        }
        
        ~ThreadData() {
            delete[] read_buffer;
        }
    };
    
    std::vector<std::unique_ptr<ThreadData>> thread_data;
    std::atomic<size_t> current_block{0};
    std::vector<std::pair<size_t, size_t>> file_blocks; // start, size
    std::atomic<size_t> total_sequences{0};
    std::atomic<size_t> sequences_processed{0};
    
    // Lookup tables
    static constexpr uint8_t char_to_bits[256] = {
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
    
    static constexpr char bits_to_char[4] = {'A', 'C', 'G', 'T'};
    
    // Оптимизированное преобразование с проверкой валидности
    inline bool string_to_kmer_fast(const char* seq, size_t pos, kmer_t& kmer) {
        kmer = 0;
        const char* p = seq + pos;
        
        // Развернутый цикл для малых k
        if (k <= 16) {
            for (size_t i = 0; i < k; ++i) {
                uint8_t bits = char_to_bits[(uint8_t)p[i]];
                if (bits == 4) return false;
                kmer = (kmer << 2) | bits;
            }
        } else {
            // Для больших k используем обычный цикл
            for (size_t i = 0; i < k; ++i) {
                uint8_t bits = char_to_bits[(uint8_t)p[i]];
                if (bits == 4) return false;
                kmer = (kmer << 2) | bits;
            }
        }
        return true;
    }
    
    // Быстрый reverse complement
    inline kmer_t reverse_complement_fast(kmer_t kmer) {
        // Используем bit manipulation tricks
        kmer = ((kmer & 0xAAAAAAAAAAAAAAAA) >> 1) | ((kmer & 0x5555555555555555) << 1);
        kmer = ((kmer & 0xCCCCCCCCCCCCCCCC) >> 2) | ((kmer & 0x3333333333333333) << 2);
        kmer = ((kmer & 0xF0F0F0F0F0F0F0F0) >> 4) | ((kmer & 0x0F0F0F0F0F0F0F0F) << 4);
        kmer = ((kmer & 0xFF00FF00FF00FF00) >> 8) | ((kmer & 0x00FF00FF00FF00FF) << 8);
        kmer = ((kmer & 0xFFFF0000FFFF0000) >> 16) | ((kmer & 0x0000FFFF0000FFFF) << 16);
        kmer = (kmer >> 32) | (kmer << 32);
        
        // Complement
        kmer = ~kmer;
        
        // Shift to correct position
        return kmer >> (64 - 2 * k);
    }
    
    inline kmer_t get_canonical(kmer_t kmer) {
        if (!use_canonical) return kmer;
        kmer_t rc = reverse_complement_fast(kmer);
        return (kmer < rc) ? kmer : rc;
    }
    
    void process_sequence_batch(ThreadData& td) {
        for (const auto& seq : td.batch) {
            if (seq.length() < k) continue;
            
            const char* s = seq.c_str();
            size_t len = seq.length();
            
            for (size_t i = 0; i <= len - k; ++i) {
                kmer_t kmer;
                if (!string_to_kmer_fast(s, i, kmer)) continue;
                
                kmer_t canonical = get_canonical(kmer);
                td.local_kmers[canonical]++;
                td.kmers_processed++;
            }
            td.sequences_processed++;
        }
        
        sequences_processed.fetch_add(td.batch.size());
        td.batch.clear();
    }
    
    void worker_thread(size_t thread_id) {
        auto& td = *thread_data[thread_id];
        
        while (true) {
            size_t block_idx = current_block.fetch_add(1);
            if (block_idx >= file_blocks.size()) break;
            
            auto [start_pos, block_size] = file_blocks[block_idx];
            
            // Читаем блок файла
            std::ifstream file("temp_input.dat", std::ios::binary);
            file.seekg(start_pos);
            file.read(td.read_buffer, block_size);
            file.close();
            
            // Парсим последовательности из блока
            std::string current_seq;
            bool in_sequence = false;
            
            for (size_t i = 0; i < block_size; ++i) {
                char c = td.read_buffer[i];
                
                if (c == '>') {
                    if (!current_seq.empty()) {
                        td.batch.push_back(std::move(current_seq));
                        current_seq.clear();
                        
                        if (td.batch.size() >= BATCH_SIZE) {
                            process_sequence_batch(td);
                        }
                    }
                    in_sequence = false;
                } else if (c == '\n') {
                    if (!in_sequence) {
                        in_sequence = true;
                    }
                } else if (in_sequence && c != '\r') {
                    current_seq += c;
                }
            }
            
            if (!current_seq.empty()) {
                td.batch.push_back(std::move(current_seq));
            }
            
            if (!td.batch.empty()) {
                process_sequence_batch(td);
            }
        }
    }
    
    void prepare_file_blocks(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary | std::ios::ate);
        size_t file_size = file.tellg();
        file.close();
        
        // Делим файл на блоки для параллельной обработки
        size_t block_size = std::max(size_t(1024 * 1024), file_size / (num_threads * 4));
        
        for (size_t pos = 0; pos < file_size; pos += block_size) {
            size_t actual_size = std::min(block_size, file_size - pos);
            file_blocks.emplace_back(pos, actual_size);
        }
    }

public:
    OptimizedKmerCounter(size_t k_value, size_t threads, size_t min_count_filter = 1, bool canonical = true)
        : k(k_value), num_threads(threads), min_count(min_count_filter), use_canonical(canonical) {
        
        thread_data.reserve(num_threads);
        for (size_t i = 0; i < num_threads; ++i) {
            thread_data.emplace_back(std::make_unique<ThreadData>());
        }
    }
    
    void count_kmers_from_file(const std::string& filename) {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Сначала копируем файл в память или создаем memory-mapped file
        std::cout << "Preparing file for parallel processing..." << std::endl;
        
        // Простой подход - читаем весь файл в память
        std::ifstream input(filename, std::ios::binary | std::ios::ate);
        size_t file_size = input.tellg();
        input.seekg(0);
        
        std::vector<char> file_buffer(file_size);
        input.read(file_buffer.data(), file_size);
        input.close();
        
        // Считаем последовательности и делим на блоки
        std::vector<size_t> sequence_starts;
        for (size_t i = 0; i < file_size; ++i) {
            if (file_buffer[i] == '>') {
                sequence_starts.push_back(i);
            }
        }
        total_sequences = sequence_starts.size();
        
        std::cout << "Found " << total_sequences << " sequences" << std::endl;
        std::cout << "Starting parallel k-mer counting with " << num_threads << " threads..." << std::endl;
        
        // Делим последовательности между потоками
        size_t seqs_per_thread = (sequence_starts.size() + num_threads - 1) / num_threads;
        
        std::vector<std::thread> workers;
        
        // Каждый поток обрабатывает свой диапазон последовательностей
        for (size_t t = 0; t < num_threads; ++t) {
            workers.emplace_back([this, t, &file_buffer, &sequence_starts, seqs_per_thread]() {
                auto& td = *thread_data[t];
                
                size_t start_idx = t * seqs_per_thread;
                size_t end_idx = std::min(start_idx + seqs_per_thread, sequence_starts.size());
                
                for (size_t i = start_idx; i < end_idx; ++i) {
                    size_t seq_start = sequence_starts[i];
                    size_t seq_end = (i + 1 < sequence_starts.size()) ? 
                                    sequence_starts[i + 1] : file_buffer.size();
                    
                    // Найти начало последовательности (после заголовка)
                    while (seq_start < seq_end && file_buffer[seq_start] != '\n') seq_start++;
                    seq_start++; // Пропускаем \n
                    
                    // Собираем последовательность
                    std::string sequence;
                    sequence.reserve(seq_end - seq_start);
                    
                    for (size_t j = seq_start; j < seq_end; ++j) {
                        char c = file_buffer[j];
                        if (c != '\n' && c != '\r' && c != '>') {
                            sequence += c;
                        }
                        if (c == '>') break;
                    }
                    
                    // Обрабатываем последовательность
                    if (sequence.length() >= k) {
                        const char* s = sequence.c_str();
                        size_t len = sequence.length();
                        
                        for (size_t pos = 0; pos <= len - k; ++pos) {
                            kmer_t kmer;
                            if (!string_to_kmer_fast(s, pos, kmer)) continue;
                            
                            kmer_t canonical = get_canonical(kmer);
                            td.local_kmers[canonical]++;
                            td.kmers_processed++;
                        }
                    }
                    
                    td.sequences_processed++;
                    sequences_processed.fetch_add(1);
                    
                    // Периодический вывод прогресса
                    if (td.sequences_processed % 10000 == 0) {
                        size_t total_proc = sequences_processed.load();
                        double progress = 100.0 * total_proc / total_sequences;
                        std::cout << "\rProgress: " << std::fixed << std::setprecision(1) 
                                 << progress << "% (" << total_proc << "/" << total_sequences << ")" 
                                 << std::flush;
                    }
                }
            });
        }
        
        // Ждем завершения всех потоков
        for (auto& t : workers) {
            t.join();
        }
        
        std::cout << "\nMerging results from all threads..." << std::endl;
        
        // Объединяем результаты
        FastMap<kmer_t, size_t> final_kmers;
        final_kmers.reserve(10000000); // Резервируем место
        
        for (const auto& td : thread_data) {
            for (const auto& [kmer, count] : td->local_kmers) {
                final_kmers[kmer] += count;
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "\nProcessing completed in " << duration.count() << " ms" << std::endl;
        std::cout << "Processing rate: " << (sequences_processed.load() * 1000.0 / duration.count()) 
                  << " sequences/second" << std::endl;
        
        // Статистика
        size_t total_kmers = 0;
        for (const auto& td : thread_data) {
            total_kmers += td->kmers_processed;
        }
        std::cout << "Total k-mers processed: " << total_kmers << std::endl;
        std::cout << "Unique k-mers found: " << final_kmers.size() << std::endl;
        
        // Сохранение результатов...
        save_results(final_kmers, "output.txt");
    }
    
    void save_results(const FastMap<kmer_t, size_t>& kmers, const std::string& output_file) {
        std::vector<std::pair<kmer_t, size_t>> sorted_kmers;
        
        for (const auto& [kmer, count] : kmers) {
            if (count >= min_count) {
                sorted_kmers.emplace_back(kmer, count);
            }
        }
        
        // Сортировка по частоте
        std::sort(sorted_kmers.begin(), sorted_kmers.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        std::ofstream out(output_file);
        for (const auto& [kmer, count] : sorted_kmers) {
            out << kmer_to_string(kmer) << "\t" << count << "\n";
        }
        out.close();
        
        std::cout << "Results saved to " << output_file << std::endl;
    }
    
    std::string kmer_to_string(kmer_t kmer) {
        std::string result(k, 'N');
        for (int i = k - 1; i >= 0; --i) {
            result[i] = bits_to_char[kmer & 3];
            kmer >>= 2;
        }
        return result;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <k> <output_file> [-t threads] [-m min_count]" << std::endl;
        return 1;
    }
    
    std::string input_file = argv[1];
    size_t k = std::stoul(argv[2]);
    std::string output_file = argv[3];
    
    size_t threads = std::thread::hardware_concurrency();
    size_t min_count = 1;
    
    // Парсинг дополнительных аргументов
    for (int i = 4; i < argc; i++) {
        if (std::string(argv[i]) == "-t" && i + 1 < argc) {
            threads = std::stoul(argv[++i]);
        } else if (std::string(argv[i]) == "-m" && i + 1 < argc) {
            min_count = std::stoul(argv[++i]);
        }
    }
    
    OptimizedKmerCounter counter(k, threads, min_count, true);
    counter.count_kmers_from_file(input_file);
    
    return 0;
}