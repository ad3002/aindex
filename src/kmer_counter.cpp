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
#include <queue>
#include <condition_variable>
#include <cctype>
#include <cstdint>
#include <bitset>

// Для k <= 32 используем uint64_t, для больших - uint128_t или вектор
using kmer_t = uint64_t;

class KmerCounter {
private:
    size_t k;
    size_t num_threads;
    size_t min_count;
    bool use_canonical;
    size_t max_k_supported = 32; // Для uint64_t максимум 32 нуклеотида
    
    std::mutex result_mutex;
    std::mutex queue_mutex;
    std::condition_variable cv;
    std::queue<std::string> sequence_queue;
    std::atomic<bool> done_reading{false};
    std::unordered_map<kmer_t, std::atomic<size_t>> kmer_counts;
    std::mutex kmer_mutex;
    
    // 2-битное кодирование: A=00, C=01, G=10, T/U=11
    static constexpr uint8_t char_to_bits[256] = {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 0-15
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 16-31
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 32-47
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 48-63
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 64-79  (@, A, B, C, D, E, F, G...)
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 80-95  (T, U)
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 96-111  (a, b, c, d, e, f, g...)
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 112-127 (t, u)
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 128-143
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 144-159
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 160-175
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 176-191
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 192-207
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 208-223
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 224-239
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4  // 240-255
    };
    
    static constexpr char bits_to_char[4] = {'A', 'C', 'G', 'T'};
    
    // Преобразование строки в 2-битное представление
    bool string_to_kmer(const std::string& seq, size_t pos, kmer_t& kmer) {
        kmer = 0;
        for (size_t i = 0; i < k; ++i) {
            uint8_t bits = char_to_bits[(uint8_t)seq[pos + i]];
            if (bits == 4) return false; // Недопустимый символ (N или другой)
            kmer = (kmer << 2) | bits;
        }
        return true;
    }
    
    // Преобразование 2-битного представления обратно в строку
    std::string kmer_to_string(kmer_t kmer) {
        std::string result(k, 'N');
        for (int i = k - 1; i >= 0; --i) {
            result[i] = bits_to_char[kmer & 3];
            kmer >>= 2;
        }
        return result;
    }
    
    // Получение обратного комплемента в 2-битном представлении
    kmer_t reverse_complement_bits(kmer_t kmer) {
        kmer_t rc = 0;
        for (size_t i = 0; i < k; ++i) {
            uint8_t bits = kmer & 3;
            // Комплемент: A(00)↔T(11), C(01)↔G(10)
            bits = 3 - bits;
            rc = (rc << 2) | bits;
            kmer >>= 2;
        }
        return rc;
    }
    
    // Получение канонической формы k-мера
    kmer_t get_canonical(kmer_t kmer) {
        if (!use_canonical) return kmer;
        kmer_t rc = reverse_complement_bits(kmer);
        return (kmer < rc) ? kmer : rc;
    }

    void process_sequence(const std::string& seq) {
        if (seq.length() < k) return;
        
        // Локальный набор для уменьшения блокировок
        std::unordered_map<kmer_t, size_t> local_kmers;
        
        for (size_t i = 0; i <= seq.length() - k; ++i) {
            kmer_t kmer;
            if (!string_to_kmer(seq, i, kmer)) continue; // Пропускаем k-меры с N
            
            kmer_t canonical_kmer = get_canonical(kmer);
            local_kmers[canonical_kmer]++;
        }
        
        // Обновляем глобальную карту
        std::lock_guard<std::mutex> lock(kmer_mutex);
        for (const auto& [kmer, count] : local_kmers) {
            kmer_counts[kmer].fetch_add(count);
        }
    }

    void worker_thread() {
        while (true) {
            std::string seq;
            
            {
                std::unique_lock<std::mutex> lock(queue_mutex);
                cv.wait(lock, [this] { return !sequence_queue.empty() || done_reading; });
                
                if (sequence_queue.empty() && done_reading) {
                    break;
                }
                
                if (!sequence_queue.empty()) {
                    seq = std::move(sequence_queue.front());
                    sequence_queue.pop();
                }
            }
            
            if (!seq.empty()) {
                process_sequence(seq);
            }
        }
    }

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
        
        // Добавляем последнюю последовательность
        if (!sequence.empty()) {
            std::lock_guard<std::mutex> lock(queue_mutex);
            sequence_queue.push(sequence);
            cv.notify_one();
        }
    }

    void read_fastq_file(std::ifstream& input) {
        std::string line;
        int line_number = 0;
        
        while (std::getline(input, line)) {
            if (line_number % 4 == 1) { // Строка с последовательностью
                if (!line.empty()) {
                    std::lock_guard<std::mutex> lock(queue_mutex);
                    sequence_queue.push(line);
                    cv.notify_one();
                }
            }
            line_number++;
        }
    }

public:
    KmerCounter(size_t k_value, size_t threads = std::thread::hardware_concurrency(), 
                size_t min_count_filter = 1, bool canonical = true) 
        : k(k_value), num_threads(threads), min_count(min_count_filter), 
          use_canonical(canonical) {
        if (num_threads == 0) num_threads = 1;
        if (k > max_k_supported) {
            throw std::runtime_error("K-mer size > " + std::to_string(max_k_supported) + 
                                   " not supported with current implementation");
        }
    }

    void count_kmers_from_file(const std::string& filename) {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Запускаем рабочие потоки
        std::vector<std::thread> workers;
        for (size_t i = 0; i < num_threads; ++i) {
            workers.emplace_back(&KmerCounter::worker_thread, this);
        }
        
        // Определяем формат файла и читаем
        std::ifstream input(filename);
        if (!input.is_open()) {
            std::cerr << "Error: Cannot open file " << filename << std::endl;
            done_reading = true;
            cv.notify_all();
            for (auto& t : workers) t.join();
            return;
        }
        
        FileFormat format = detect_format(filename);
        input.close();
        input.open(filename);
        
        std::cout << "File format detected: ";
        switch (format) {
            case FileFormat::FASTA: std::cout << "FASTA" << std::endl; break;
            case FileFormat::FASTQ: std::cout << "FASTQ" << std::endl; break;
            case FileFormat::PLAIN: std::cout << "Plain text" << std::endl; break;
        }
        
        std::cout << "Using 2-bit encoding (memory usage reduced ~4x)" << std::endl;
        
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
        done_reading = true;
        cv.notify_all();
        
        // Ждем завершения всех потоков
        for (auto& t : workers) {
            t.join();
        }
        
        // Фильтруем по минимальной частоте
        if (min_count > 1) {
            auto it = kmer_counts.begin();
            while (it != kmer_counts.end()) {
                if (it->second.load() < min_count) {
                    it = kmer_counts.erase(it);
                } else {
                    ++it;
                }
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Processing completed in " << duration.count() << " ms" << std::endl;
        std::cout << "Using " << (use_canonical ? "canonical" : "non-canonical") << " k-mers" << std::endl;
        if (min_count > 1) {
            std::cout << "Applied minimum count filter: " << min_count << std::endl;
        }
        std::cout << "Found " << kmer_counts.size() << " unique " << k << "-mers" << std::endl;
        
        // Оценка использования памяти
        size_t memory_mb = (kmer_counts.size() * (sizeof(kmer_t) + sizeof(std::atomic<size_t>) + 32)) / (1024 * 1024);
        std::cout << "Estimated memory usage: ~" << memory_mb << " MB" << std::endl;
    }

    void save_kmers(const std::string& output_file, bool with_counts = true) {
        std::ofstream out(output_file);
        if (!out.is_open()) {
            std::cerr << "Error: Cannot create output file " << output_file << std::endl;
            return;
        }
        
        // Собираем k-меры в вектор для сортировки
        std::vector<std::pair<kmer_t, size_t>> kmers;
        for (const auto& [kmer, count] : kmer_counts) {
            size_t c = count.load();
            if (c >= min_count) {
                kmers.emplace_back(kmer, c);
            }
        }
        
        // Сортировка по частоте (убывание)
        std::sort(kmers.begin(), kmers.end(), 
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        // Запись результатов
        for (const auto& [kmer, count] : kmers) {
            std::string kmer_str = kmer_to_string(kmer);
            if (with_counts) {
                out << kmer_str << "\t" << count << "\n";
            } else {
                out << kmer_str << "\n";
            }
        }
        
        out.close();
        std::cout << "Results saved to " << output_file << std::endl;
    }

    void save_kmers_binary(const std::string& output_file) {
        std::ofstream out(output_file, std::ios::binary);
        if (!out.is_open()) {
            std::cerr << "Error: Cannot create binary output file " << output_file << std::endl;
            return;
        }
        
        // Считаем количество k-меров после фильтрации
        size_t num_kmers = 0;
        for (const auto& [kmer, count] : kmer_counts) {
            if (count.load() >= min_count) num_kmers++;
        }
        
        // Записываем заголовок
        out.write(reinterpret_cast<const char*>(&num_kmers), sizeof(size_t));
        out.write(reinterpret_cast<const char*>(&k), sizeof(size_t));
        
        // Записываем k-меры в компактном формате
        for (const auto& [kmer, count] : kmer_counts) {
            size_t c = count.load();
            if (c >= min_count) {
                out.write(reinterpret_cast<const char*>(&kmer), sizeof(kmer_t));
                out.write(reinterpret_cast<const char*>(&c), sizeof(size_t));
            }
        }
        
        out.close();
        std::cout << "Binary results saved to " << output_file 
                  << " (2-bit encoded, ~" << (num_kmers * (sizeof(kmer_t) + sizeof(size_t))) / 1024 / 1024 
                  << " MB)" << std::endl;
    }

    // Сохранение в компактном бинарном формате для perfect hash
    void save_compact_binary(const std::string& output_file) {
        std::ofstream out(output_file, std::ios::binary);
        if (!out.is_open()) {
            std::cerr << "Error: Cannot create compact binary file " << output_file << std::endl;
            return;
        }
        
        // Заголовок: магическое число, версия, k, количество k-меров
        uint32_t magic = 0x4B4D4552; // "KMER"
        uint32_t version = 1;
        uint32_t k32 = k;
        uint64_t num_kmers = 0;
        
        for (const auto& [kmer, count] : kmer_counts) {
            if (count.load() >= min_count) num_kmers++;
        }
        
        out.write(reinterpret_cast<const char*>(&magic), sizeof(uint32_t));
        out.write(reinterpret_cast<const char*>(&version), sizeof(uint32_t));
        out.write(reinterpret_cast<const char*>(&k32), sizeof(uint32_t));
        out.write(reinterpret_cast<const char*>(&num_kmers), sizeof(uint64_t));
        
        // Сортируем k-меры по значению для лучшей локальности
        std::vector<std::pair<kmer_t, uint32_t>> sorted_kmers;
        for (const auto& [kmer, count] : kmer_counts) {
            size_t c = count.load();
            if (c >= min_count) {
                sorted_kmers.emplace_back(kmer, std::min(c, (size_t)UINT32_MAX));
            }
        }
        std::sort(sorted_kmers.begin(), sorted_kmers.end());
        
        // Записываем k-меры
        for (const auto& [kmer, count] : sorted_kmers) {
            out.write(reinterpret_cast<const char*>(&kmer), sizeof(kmer_t));
            if (version == 1) { // В версии 1 сохраняем счетчики
                out.write(reinterpret_cast<const char*>(&count), sizeof(uint32_t));
            }
        }
        
        out.close();
        std::cout << "Compact binary saved to " << output_file << ".cbin" << std::endl;
    }

    // Метод для получения статистики
    void print_statistics() {
        if (kmer_counts.empty()) {
            std::cout << "No k-mers found" << std::endl;
            return;
        }
        
        size_t total_kmers = 0;
        size_t max_count = 0;
        size_t singleton_count = 0;
        size_t filtered_count = 0;
        
        std::vector<size_t> counts;
        for (const auto& [kmer, count] : kmer_counts) {
            size_t c = count.load();
            if (c >= min_count) {
                counts.push_back(c);
                total_kmers += c;
                max_count = std::max(max_count, c);
                if (c == 1) singleton_count++;
            } else {
                filtered_count++;
            }
        }
        
        // Медиана
        size_t median = 0;
        if (!counts.empty()) {
            std::sort(counts.begin(), counts.end());
            median = counts[counts.size() / 2];
        }
        
        std::cout << "\n=== K-mer Statistics ===" << std::endl;
        std::cout << "Total k-mers: " << total_kmers << std::endl;
        std::cout << "Unique k-mers: " << kmer_counts.size() - filtered_count << std::endl;
        if (filtered_count > 0) {
            std::cout << "Filtered k-mers: " << filtered_count << std::endl;
        }
        std::cout << "Singleton k-mers: " << singleton_count << std::endl;
        std::cout << "Max k-mer frequency: " << max_count << std::endl;
        std::cout << "Median frequency: " << median << std::endl;
        if (!counts.empty()) {
            std::cout << "Average frequency: " << (double)total_kmers / counts.size() << std::endl;
        }
        
        // Теоретический максимум k-меров
        size_t theoretical_max = 1ULL << (2 * k);
        if (theoretical_max > 0) {
            double coverage = 100.0 * (kmer_counts.size() - filtered_count) / theoretical_max;
            std::cout << "K-mer space coverage: " << coverage << "%" << std::endl;
        }
    }
};

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <k> <output_file> [options]" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "  -t <threads>     Number of threads (default: auto)" << std::endl;
        std::cerr << "  -m <min_count>   Minimum k-mer count (default: 1)" << std::endl;
        std::cerr << "  -c               Use canonical k-mers (default: yes)" << std::endl;
        std::cerr << "  -n               Don't use canonical k-mers" << std::endl;
        std::cerr << "  -b               Save compact binary format (.cbin)" << std::endl;
        std::cerr << "\nNote: Maximum k-mer size is 32 for 64-bit encoding" << std::endl;
        std::cerr << "\nExample: " << argv[0] << " genome.fasta 31 kmers.txt -t 8 -m 2" << std::endl;
        return 1;
    }
    
    std::string input_file = argv[1];
    size_t k = std::stoul(argv[2]);
    std::string output_file = argv[3];
    
    // Параметры по умолчанию
    size_t threads = std::thread::hardware_concurrency();
    size_t min_count = 1;
    bool use_canonical = true;
    bool save_compact = false;
    
    // Парсинг опций
    for (int i = 4; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-t" && i + 1 < argc) {
            threads = std::stoul(argv[++i]);
        } else if (arg == "-m" && i + 1 < argc) {
            min_count = std::stoul(argv[++i]);
        } else if (arg == "-c") {
            use_canonical = true;
        } else if (arg == "-n") {
            use_canonical = false;
        } else if (arg == "-b") {
            save_compact = true;
        }
    }
    
    try {
        std::cout << "K-mer Counter Configuration:" << std::endl;
        std::cout << "Input file: " << input_file << std::endl;
        std::cout << "K-mer size: " << k << std::endl;
        std::cout << "Output file: " << output_file << std::endl;
        std::cout << "Threads: " << threads << std::endl;
        std::cout << "Min count filter: " << min_count << std::endl;
        std::cout << "Canonical k-mers: " << (use_canonical ? "yes" : "no") << std::endl;
        std::cout << "2-bit encoding: enabled (4x memory reduction)" << std::endl << std::endl;
        
        KmerCounter counter(k, threads, min_count, use_canonical);
        counter.count_kmers_from_file(input_file);
        counter.print_statistics();
        
        // Сохраняем в текстовом формате
        counter.save_kmers(output_file, true);
        
        // Сохраняем в бинарном формате
        std::string binary_output = output_file + ".bin";
        counter.save_kmers_binary(binary_output);
        
        // Опционально сохраняем компактный бинарный формат
        if (save_compact) {
            std::string compact_output = output_file + ".cbin";
            counter.save_compact_binary(compact_output);
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}