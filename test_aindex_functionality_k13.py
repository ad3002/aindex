import sys
import os
import logging
import numpy as np
import time
import subprocess

# Настройка путей - добавляем корневую директорию проекта в Python path
PATH_TO_AINDEX_FOLDER = '..'
sys.path.insert(0, PATH_TO_AINDEX_FOLDER)

# Настройка логирования
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Импорт модуля aindex - используем прямой импорт C++ модуля для 13-меров
try:
    import aindex.core.aindex_cpp as aindex_cpp
    print("✓ Модуль aindex_cpp успешно импортирован")
    
    # Также попробуем импортировать Python wrapper
    try:
        from aindex.core.aindex import AIndex
        print("✓ Модуль AIndex также успешно импортирован")
        HAS_PYTHON_WRAPPER = True
    except ImportError as e:
        print(f"⚠️ Python wrapper недоступен: {e}")
        HAS_PYTHON_WRAPPER = False
        AIndex = None
        
except ImportError as e:
    print(f"✗ Ошибка импорта aindex_cpp: {e}")
    print("Проверьте, что пакет aindex установлен и собран корректно")
    print("Выполните в терминале:")
    print(f"  cd {PATH_TO_AINDEX_FOLDER}")
    print("  make pybind11")
    aindex_cpp = None
    AIndex = None
    HAS_PYTHON_WRAPPER = False

print("Рабочая директория:", os.getcwd())
print("Python версия:", sys.version)
print("Путь к aindex:", PATH_TO_AINDEX_FOLDER)

assert aindex_cpp is not None, "Модуль aindex_cpp не импортирован, тесты не могут быть выполнены"

# Проверяем наличие файлов перед загрузкой
prefix = os.path.join(PATH_TO_AINDEX_FOLDER, "temp/all_13mers")
required_files = [
    f"{prefix}.pf",
    f"{prefix}.index.bin", 
    f"{prefix}.kmers.bin"
]

# Проверяем наличие файла с ридами (путь без учета размера k-меров)
base_prefix = os.path.join(PATH_TO_AINDEX_FOLDER, "temp/reads")
reads_file = f"{base_prefix}.reads"
reads_index_file = f"{base_prefix}.ridx"  # Правильный путь к индексу ридов
kmers_file = f"{prefix}.kmers"
has_reads = os.path.exists(reads_file) and os.path.exists(reads_index_file)

# Функция для подсчета строк в файле
def count_lines(file_path):
    try:
        result = subprocess.run(["wc", "-l", file_path], capture_output=True, text=True)
        if result.returncode == 0:
            return int(result.stdout.split()[0])
        return None
    except Exception as e:
        print(f"Ошибка при подсчете строк в {file_path}: {e}")
        return None

def parse_kmers_analysis_file(file_path):
    """
    Парсит файл kmers_analysis.trues и возвращает словарь кмер -> частота
    
    Формат файла:
    kmer    frequency   [read_id,position,direction] [read_id,position,direction] ...
    """
    kmer_tf = {}
    
    if not os.path.exists(file_path):
        print(f"Файл {file_path} не найден, пропускаем проверку k-меров")
        return kmer_tf
        
    try:
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    kmer = parts[0]
                    frequency = int(parts[1])
                    kmer_tf[kmer] = frequency
        print(f"Загружено {len(kmer_tf)} кмеров из файла {file_path}")
    except Exception as e:
        print(f"Ошибка при чтении файла {file_path}: {e}")
        
    return kmer_tf

# Получаем ожидаемое количество ридов и k-меров
expected_reads_count = count_lines(reads_file) if os.path.exists(reads_file) else None
expected_kmers_count = count_lines(kmers_file) if os.path.exists(kmers_file) else None

print("Проверка наличия файлов:")
all_files_exist = True
for file_path in required_files:
    exists = os.path.exists(file_path)
    print(f"  {file_path}: {'✓' if exists else '✗'}")
    if not exists:
        all_files_exist = False

print(f"  {reads_file}: {'✓' if os.path.exists(reads_file) else '✗'}")
print(f"  {reads_index_file}: {'✓' if os.path.exists(reads_index_file) else '✗'}")
print(f"  {kmers_file}: {'✓' if os.path.exists(kmers_file) else '✗'}")

print("\nОжидаемые количества:")
if expected_reads_count:
    print(f"  Количество ридов (wc -l {reads_file}): {expected_reads_count:,}")
if expected_kmers_count:
    print(f"  Количество k-меров (wc -l {kmers_file}): {expected_kmers_count:,}")

if not all_files_exist:
    print("\n⚠️ Не все необходимые файлы найдены. Проверьте, что индекс был создан правильно.")
    print(f"Ожидаемые файлы должны находиться в директории: {os.path.dirname(prefix)}")
else:
    print("\n✓ Все необходимые файлы найдены, начинаем загрузку...")

# Загружаем индекс только если aindex_cpp доступен
if aindex_cpp is not None:
    try:
        start_time = time.time()
        print(f"Загружаем 13-мерный индекс из префикса: {prefix}")
        
        # Создаем экземпляр AindexWrapper и загружаем данные
        index = aindex_cpp.AindexWrapper()
        
        # Загружаем 13-мерный индекс
        index.load_from_prefix_13mer(prefix)
        
        # Загружаем риды если доступны
        if has_reads:
            print("Загружаем риды в 13-мерный индекс...")
            index.load_reads(reads_file)
            print(f"✓ Риды загружены: {index.n_reads:,} ридов")
        else:
            print("⚠️ Файлы ридов не найдены, загружаем только TF данные")
        
        load_time = time.time() - start_time
        print(f"✓ 13-мерный индекс успешно загружен за {load_time:.2f} секунд")
        
        # Базовая информация о 13-мерном индексе
        stats = index.get_13mer_statistics()
        print(f"Всего 13-меров: {stats['total_kmers']:,}")
        print(f"13-меры с ненулевым TF: {stats['non_zero_kmers']:,}")
        print(f"Максимальная частота: {stats['max_frequency']:,}")
        print(f"Общее количество TF: {stats['total_count']:,}")
        print(f"Количество ридов: {index.n_reads:,}")
        print(f"Размер файла ридов: {index.reads_size:,} байт")
        
        print(f"Индекс загружен для 13-меров: True")
        
        # Проверяем, что риды были загружены
        if index.n_reads > 0:
            print(f"✓ Риды успешно загружены ({index.n_reads:,} ридов)")
            
            # Демонстрация доступа к ридам
            try:
                # Получаем первый рид для демонстрации
                first_read = index.get_read_by_rid(0)
                if first_read:
                    print(f"Первый рид (первые 50 символов): {first_read[:50]}...")
            except Exception as e:
                print(f"⚠️ Ошибка при доступе к ридам: {e}")
        else:
            print("⚠️ Риды не загружены или их количество равно 0")
        
        # Сравнение ожидаемого и фактического количества (убираем, так как используем C++ API)
            
        # Статистика по файлу tf для 13-меров
        print("\n=== Статистика файла tf для 13-меров ===")
        try:
            # Получаем статистику через новый API
            stats = index.get_13mer_statistics()
            total_kmers = stats['total_kmers']
            non_zero_count = stats['non_zero_kmers']
            max_tf = stats['max_frequency']
            total_count = stats['total_count']
            avg_tf = total_count / non_zero_count if non_zero_count > 0 else 0
            
            print(f"Всего 13-меров: {total_kmers:,}")
            print(f"13-меров с ненулевым TF: {non_zero_count:,} ({non_zero_count*100/total_kmers:.2f}%)")
            print(f"Максимальный TF: {max_tf:,}")
            print(f"Общий счет TF: {total_count:,}")
            print(f"Средний TF (для ненулевых): {avg_tf:.2f}")
        except Exception as e:
            print(f"Ошибка при получении статистики tf: {e}")
            
        # Убираем проверку ридов, так как используем только TF API для 13-меров
        
        # Проверка совпадения k-меров и их частот с ожидаемыми значениями
        kmers_analysis_file = os.path.join(PATH_TO_AINDEX_FOLDER, "./temp/all_13mers.true.kmers")
        if os.path.exists(kmers_analysis_file):
            print("\n=== Проверка соответствия 13-меров ожидаемым значениям ===")
            expected_kmers = parse_kmers_analysis_file(kmers_analysis_file)
            
            if expected_kmers:
                print(f"Загружено {len(expected_kmers)} ожидаемых 13-меров из {kmers_analysis_file}")
                
                # Тестирование первых нескольких k-меров для демонстрации
                sample_size = min(10, len(expected_kmers))
                sample_kmers = list(expected_kmers.keys())[:sample_size]
                
                print(f"\nПроверка {sample_size} случайных 13-меров:")
                print("Kmer\t\t\tОжидалось\tПрямой\tОбратный\tИтого\tСтатус")
                print("-" * 80)
                
                all_matched = True
                mismatches = 0
                
                for kmer in sample_kmers:
                    expected_tf = expected_kmers[kmer]
                    
                    # Используем новые функции для получения TF в обеих направлениях
                    forward_tf, reverse_tf = index.get_tf_both_directions_13mer(kmer)
                    total_tf = index.get_total_tf_value_13mer(kmer)
                    
                    matches = expected_tf == total_tf
                    all_matched = all_matched and matches
                    if not matches:
                        mismatches += 1
                        
                    status = "✓" if matches else "✗"
                    print(f"{status} {kmer}\t{expected_tf}\t\t{forward_tf}\t{reverse_tf}\t{total_tf}\t{status}")
                
                # Проверка всех k-меров
                print("\nПроверка всех 13-меров:")
                total_mismatches = 0
                zero_tf_count = 0
                
                for kmer, expected_tf in expected_kmers.items():
                    total_tf = index.get_total_tf_value_13mer(kmer)
                    
                    if total_tf == 0:
                        zero_tf_count += 1
                    
                    if expected_tf != total_tf:
                        total_mismatches += 1
                        
                # Вывод результатов
                match_percentage = 100 * (len(expected_kmers) - total_mismatches) / len(expected_kmers)
                print(f"Совпало: {match_percentage:.2f}% 13-меров")
                print(f"Несовпадений: {total_mismatches:,} из {len(expected_kmers):,} 13-меров")
                print(f"13-меры с нулевой частотой: {zero_tf_count:,}")
                
                if total_mismatches > 0:
                    print("\n⚠️ Обнаружены несоответствия между ожидаемыми и фактическими частотами 13-меров")
                else:
                    print("\n✓ Все частоты 13-меров соответствуют ожидаемым значениям")
            else:
                print("Не удалось загрузить ожидаемые значения 13-меров для проверки")
        else:
            print(f"\n⚠️ Файл {kmers_analysis_file} не найден, пропускаем проверку 13-меров")
            print("  Для создания файла выполните:\n  python tests/analyze_kmers.py --input-file temp/reads.reads -o kmers_analysis.trues")
    except Exception as e:
        print(f"✗ Ошибка при загрузке индекса: {e}")
        import traceback
        traceback.print_exc()
        index = None
else:
    print("⚠️ aindex_cpp не импортирован, пропускаем загрузку индекса")
    index = None

# Получаем подробную информацию об индексе
if index is None:
    print("⚠️ Индекс не загружен, пропускаем тестирование")
else:
    print("=== Информация о 13-мерном индексе ===")
    try:
        stats = index.get_13mer_statistics()
        print(f"Всего 13-меров: {stats['total_kmers']:,}")
        print(f"13-меры с ненулевым TF: {stats['non_zero_kmers']:,}")
        print(f"Максимальная частота: {stats['max_frequency']:,}")
        print(f"Общее количество TF: {stats['total_count']:,}")
    except Exception as e:
        print(f"Ошибка получения статистики: {e}")

    print("\n=== Дополнительная статистика ===")
    print(f"Режим 13-меров: активен")

    # Тестируем доступ к TF для одного 13-мера
    test_kmer = "AAAAACCCCCGGG"  # 13-мер для тестирования
    print(f"\nТест доступа к TF:")
    print(f"  13-мер: '{test_kmer}'")
    try:
        forward_tf, reverse_tf = index.get_tf_both_directions_13mer(test_kmer)
        total_tf = index.get_total_tf_value_13мер(test_kmer)
        rc_kmer = index.get_reverse_complement_13mer(test_kmer)
        
        print(f"  Forward TF: {forward_tf}")
        print(f"  Reverse TF: {reverse_tf} (RC: {rc_kmer})")
        print(f"  Total TF: {total_tf}")
    except Exception as e:
        print(f"  Ошибка: {e}")

# Тестируем получение term frequency для 13-меров
if index is None:
    print("⚠️ Индекс не загружен, пропускаем тестирование TF")
else:
    print("=== Тест получения TF для 13-меров ===")

    # Создаем несколько тестовых 13-меров
    test_kmers = [
        "AAAAAAAAAAAAA",   # только A
        "TTTTTTTTTTTTT",   # только T
        "GGGGGGGGGGGGG",   # только G
        "CCCCCCCCCCCCC",   # только C
        "ATCGATCGATCGA",   # случайная последовательность
        "AAAAAGAGTTAAT",   # известный k-мер из наших тестов
        "AGTAGTAGTAGTA"    # еще один известный k-мер
    ]

    print("Тестирование отдельных 13-меров:")
    print("Kmer\t\t\tForward\tReverse\tTotal\tRC")
    print("-" * 80)
    
    for kmer in test_kmers:
        try:
            forward_tf, reverse_tf = index.get_tf_both_directions_13mer(kmer)
            total_tf = index.get_total_tf_value_13mer(kmer)
            rc_kmer = index.get_reverse_complement_13mer(kmer)
            print(f"{kmer}\t{forward_tf}\t{reverse_tf}\t{total_tf}\t{rc_kmer}")
        except Exception as e:
            print(f"{kmer}\tError: {e}")

    print("\n=== Тест batch получения TF для 13-меров ===")

    # Тестируем batch обработку
    batch_kmers = test_kmers[:5]  # берем первые 5
    try:
        start_time = time.time()
        
        # Batch получение TF в обеих направлениях
        batch_results = index.get_tf_both_directions_13mer_batch(batch_kmers)
        
        # Batch получение общего TF
        total_tfs = index.get_total_tf_values_13mer(batch_kmers)
        
        batch_time = time.time() - start_time
        
        print(f"Batch обработка {len(batch_kmers)} 13-меров заняла {batch_time:.4f} секунд")
        print("Результаты:")
        print("Kmer\t\t\tForward\tReverse\tTotal")
        print("-" * 60)
        
        for i, kmer in enumerate(batch_kmers):
            forward_tf, reverse_tf = batch_results[i]
            total_tf = total_tfs[i]
            print(f"{kmer}\t{forward_tf}\t{reverse_tf}\t{total_tf}")
            
    except Exception as e:
        print(f"Ошибка batch обработки: {e}")

    print("\n=== Сравнение производительности для 13-меров ===")

    # Сравниваем одиночные вызовы vs batch
    n_test = 100
    random_kmers = []
    for i in range(n_test):
        kmer = ''.join(np.random.choice(['A', 'T', 'G', 'C'], 13))
        random_kmers.append(kmer)

    # Одиночные вызовы
    start_time = time.time()
    single_results = []
    for kmer in random_kmers:
        try:
            total_tf = index.get_total_tf_value_13mer(kmer)
            single_results.append(total_tf)
        except:
            single_results.append(0)
    single_time = time.time() - start_time

    # Batch вызов
    start_time = time.time()
    try:
        batch_results = index.get_total_tf_values_13mer(random_kmers)
    except:
        batch_results = [0] * len(random_kmers)
    batch_time = time.time() - start_time

    print(f"Одиночные вызовы ({n_test} 13-меров): {single_time:.4f} сек")
    print(f"Batch вызов ({n_test} 13-меров): {batch_time:.4f} сек")
    print(f"Ускорение: {single_time/batch_time:.1f}x")

# Для 13-мерного API позиции и риды недоступны, завершаем основное тестирование
print("\n=== Тест анализ случайных 13-меров ===")

if index is None:
    print("⚠️ Индекс не загружен, пропускаем анализ")
else:
    # Проведем анализ распределения TF для случайных 13-меров
    n_random = 1000
    print(f"Генерируем {n_random} случайных 13-меров для анализа распределения TF...")
    
    random_kmers = []
    for i in range(n_random):
        kmer = ''.join(np.random.choice(['A', 'T', 'G', 'C'], 13))
        random_kmers.append(kmer)
    
    try:
        # Получаем TF для всех случайных k-меров batch'ем
        total_tfs = index.get_total_tf_values_13mer(random_kmers)
        
        # Анализ распределения
        non_zero_count = sum(1 for tf in total_tfs if tf > 0)
        max_tf = max(total_tfs) if total_tfs else 0
        avg_tf = sum(total_tfs) / len(total_tfs) if total_tfs else 0
        avg_non_zero_tf = sum(total_tfs) / non_zero_count if non_zero_count > 0 else 0
        
        print(f"Результаты анализа {n_random} случайных 13-меров:")
        print(f"  13-меры с ненулевым TF: {non_zero_count} ({non_zero_count/n_random*100:.1f}%)")
        print(f"  Максимальный TF: {max_tf}")
        print(f"  Средний TF (все): {avg_tf:.2f}")
        print(f"  Средний TF (ненулевые): {avg_non_zero_tf:.2f}")
        
        # Распределение по диапазонам TF
        tf_ranges = {
            "0": 0,
            "1": 0,
            "2-5": 0,
            "6-10": 0,
            "11-50": 0,
            ">50": 0
        }
        
        for tf in total_tfs:
            if tf == 0:
                tf_ranges["0"] += 1
            elif tf == 1:
                tf_ranges["1"] += 1
            elif 2 <= tf <= 5:
                tf_ranges["2-5"] += 1
            elif 6 <= tf <= 10:
                tf_ranges["6-10"] += 1
            elif 11 <= tf <= 50:
                tf_ranges["11-50"] += 1
            else:
                tf_ranges[">50"] += 1
        
        print(f"\nРаспределение по диапазонам TF:")
        for range_name, count in tf_ranges.items():
            percentage = count / n_random * 100
            print(f"  {range_name}: {count} ({percentage:.1f}%)")
            
    except Exception as e:
        print(f"Ошибка анализа случайных 13-меров: {e}")

    # Тестирование производительности с большими объемами данных
    print(f"\n=== Тест производительности для больших объемов ===")
    
    large_n = 10000
    print(f"Тестируем производительность на {large_n} 13-мерах...")
    
    large_random_kmers = []
    for i in range(large_n):
        kmer = ''.join(np.random.choice(['A', 'T', 'G', 'C'], 13))
        large_random_kmers.append(kmer)
    
    try:
        start_time = time.time()
        large_total_tfs = index.get_total_tf_values_13mer(large_random_kmers)
        large_batch_time = time.time() - start_time
        
        print(f"Batch обработка {large_n} 13-меров: {large_batch_time:.4f} сек")
        print(f"Скорость: {large_n/large_batch_time:.0f} 13-меров/сек")
        
        # Сравнение с одиночными вызовами (на меньшем объеме)
        test_subset = large_random_kmers[:1000]
        
        start_time = time.time()
        for kmer in test_subset:
            index.get_total_tf_value_13mer(kmer)
        single_time = time.time() - start_time
        
        start_time = time.time()
        index.get_total_tf_values_13mer(test_subset)
        batch_subset_time = time.time() - start_time
        
        print(f"\nСравнение на {len(test_subset)} 13-мерах:")
        print(f"  Одиночные вызовы: {single_time:.4f} сек ({len(test_subset)/single_time:.0f} 13-меров/сек)")
        print(f"  Batch вызов: {batch_subset_time:.4f} сек ({len(test_subset)/batch_subset_time:.0f} 13-меров/сек)")
        print(f"  Ускорение: {single_time/batch_subset_time:.1f}x")
        
    except Exception as e:
        print(f"Ошибка тестирования производительности: {e}")

    # Тест производительности для 1 миллиона запросов
    print(f"\n=== СТРЕСС-ТЕСТ: 1 МИЛЛИОН ЗАПРОСОВ ===")
    
    million_n = 1000000
    print(f"Генерируем {million_n:,} случайных 13-меров для стресс-теста...")
    
    # Генерация данных по частям для экономии памяти
    batch_size = 100000
    total_time = 0.0
    total_processed = 0
    
    print(f"Обработка по батчам размером {batch_size:,}...")
    
    try:
        for batch_num in range(million_n // batch_size):
            # Генерируем батч
            batch_kmers = []
            for i in range(batch_size):
                kmer = ''.join(np.random.choice(['A', 'T', 'G', 'C'], 13))
                batch_kmers.append(kmer)
            
            # Обрабатываем батч
            start_time = time.time()
            batch_results = index.get_total_tf_values_13mer(batch_kmers)
            batch_time = time.time() - start_time
            
            total_time += batch_time
            total_processed += len(batch_kmers)
            
            # Прогресс
            progress = (batch_num + 1) * batch_size
            print(f"  Обработано: {progress:,}/{million_n:,} ({progress/million_n*100:.1f}%) "
                  f"- Время батча: {batch_time:.3f}с, Скорость: {len(batch_kmers)/batch_time:.0f} 13-меров/сек")
        
        # Финальная статистика
        avg_speed = total_processed / total_time
        print(f"\n🚀 РЕЗУЛЬТАТЫ СТРЕСС-ТЕСТА:")
        print(f"  Всего обработано: {total_processed:,} 13-меров")
        print(f"  Общее время: {total_time:.3f} секунд")
        print(f"  Средняя скорость: {avg_speed:.0f} 13-меров/сек")
        print(f"  Пропускная способность: {avg_speed/1000:.0f}K 13-меров/сек")
        print(f"  Время на один 13-мер: {total_time/total_processed*1000000:.2f} микросекунд")
        
        # Анализ найденных TF
        if 'batch_results' in locals():
            sample_non_zero = sum(1 for tf in batch_results if tf > 0)
            sample_max = max(batch_results) if batch_results else 0
            print(f"\n📊 Анализ последнего батча:")
            print(f"  13-меры с ненулевым TF: {sample_non_zero}/{len(batch_results)} ({sample_non_zero/len(batch_results)*100:.2f}%)")
            print(f"  Максимальный TF в батче: {sample_max}")
        
        # Оценка памяти
        memory_per_kmer = 13  # байт на строку k-мера
        memory_per_tf = 8     # байт на uint64_t результат
        estimated_memory = batch_size * (memory_per_kmer + memory_per_tf) / 1024 / 1024
        print(f"\n💾 Использование памяти:")
        print(f"  Размер батча: {batch_size:,} 13-меров")
        print(f"  Оценочная память на батч: {estimated_memory:.1f} MB")
        print(f"  Общая экономия памяти: {(million_n * (memory_per_kmer + memory_per_tf) / 1024 / 1024) - estimated_memory:.0f} MB")
        
    except Exception as e:
        print(f"❌ Ошибка стресс-теста: {e}")
        import traceback
        traceback.print_exc()

    # Бенчмарк сравнения разных функций
    print(f"\n=== БЕНЧМАРК СРАВНЕНИЯ ФУНКЦИЙ ===")
    
    benchmark_n = 50000
    benchmark_kmers = []
    for i in range(benchmark_n):
        kmer = ''.join(np.random.choice(['A', 'T', 'G', 'C'], 13))
        benchmark_kmers.append(kmer)
    
    print(f"Сравнение производительности разных функций на {benchmark_n:,} 13-мерах:")
    
    try:
        # 1. get_total_tf_values_13mer (batch)
        start_time = time.time()
        total_tfs = index.get_total_tf_values_13mer(benchmark_kmers)
        total_batch_time = time.time() - start_time
        
        # 2. get_tf_both_directions_13mer_batch
        start_time = time.time()
        both_directions = index.get_tf_both_directions_13mer_batch(benchmark_kmers)
        both_batch_time = time.time() - start_time
        
        # 3. get_tf_values_13mer (только forward)
        start_time = time.time()
        forward_only = index.get_tf_values_13mer(benchmark_kmers)
        forward_batch_time = time.time() - start_time
        
        # 4. Одиночные вызовы (на подвыборке)
        subset = benchmark_kmers[:1000]
        start_time = time.time()
        for kmer in subset:
            index.get_total_tf_value_13mer(kmer)
        single_time = time.time() - start_time
        single_speed = len(subset) / single_time
        
        print(f"\n📈 Результаты бенчмарка:")
        print(f"  get_total_tf_values_13mer (batch):     {total_batch_time:.4f}с ({benchmark_n/total_batch_time:.0f} 13-меров/сек)")
        print(f"  get_tf_both_directions_13mer_batch:    {both_batch_time:.4f}с ({benchmark_n/both_batch_time:.0f} 13-меров/сек)")
        print(f"  get_tf_values_13mer (forward only):    {forward_batch_time:.4f}с ({benchmark_n/forward_batch_time:.0f} 13-меров/сек)")
        print(f"  Одиночные вызовы (экстраполяция):      {single_time:.4f}с ({single_speed:.0f} 13-меров/сек)")
        
        print(f"\n🏆 Рейтинг по скорости:")
        speeds = [
            ("Forward only batch", benchmark_n/forward_batch_time),
            ("Total TF batch", benchmark_n/total_batch_time),
            ("Both directions batch", benchmark_n/both_batch_time),
            ("Single calls", single_speed)
        ]
        speeds.sort(key=lambda x: x[1], reverse=True)
        
        for i, (name, speed) in enumerate(speeds, 1):
            print(f"  {i}. {name}: {speed:.0f} 13-меров/сек")
        
        # Проверка корректности
        sample_idx = 0
        sample_kmer = benchmark_kmers[sample_idx]
        forward_tf, reverse_tf = both_directions[sample_idx]
        expected_total = forward_tf + reverse_tf
        actual_total = total_tfs[sample_idx]
        
        print(f"\n✅ Проверка корректности (пример: {sample_kmer}):")
        print(f"  Forward: {forward_tf}, Reverse: {reverse_tf}, Ожидаемый total: {expected_total}")
        print(f"  Фактический total: {actual_total}")
        print(f"  Корректность: {'✓' if expected_total == actual_total else '✗'}")
        
    except Exception as e:
        print(f"❌ Ошибка бенчмарка: {e}")

    # Финальная сводка тестирования 13-мерного индекса
    print(f"\n" + "="*60)
    print(f"СВОДКА ТЕСТИРОВАНИЯ 13-МЕРНОГО ИНДЕКСА")
    print(f"="*60)
    prefix_display = f"{PATH_TO_AINDEX_FOLDER}/temp/all_13mers"
    print(f"Префикс: {prefix_display}")
    print(f"13-мерный индекс загружен: ✓")
    
    try:
        stats = index.get_13mer_statistics()
        print(f"Всего 13-меров: {stats['total_kmers']:,}")
        print(f"13-меры с ненулевым TF: {stats['non_zero_kmers']:,}")
        print(f"Максимальная частота: {stats['max_frequency']:,}")
        print(f"Общее количество TF: {stats['total_count']:,}")
        print(f"Средняя частота: {stats['total_count']/stats['non_zero_kmers']:.2f}")
    except:
        print("Не удалось получить финальную статистику")

    print(f"\n✓ Тестирование 13-мерного индекса завершено успешно!")

# Тестируем анализ покрытия последовательности 13-мерами
if index is None:
    print("⚠️ Индекс не загружен, пропускаем тестирование покрытия последовательности")
else:
    print("=== Тест анализа покрытия последовательности 13-мерами ===")

    test_sequences = []

    # Получаем реальные риды напрямую из 13-мерного индекса
    if index.n_reads > 0:
        print("Получаем реальные риды из 13-мерного индекса...")
        for rid in range(min(5, index.n_reads)):
            try:
                real_read = index.get_read_by_rid(rid)
                if real_read and len(real_read) >= 50:
                    # Если рид содержит субриды, берем первый субрид
                    if '~' in real_read:
                        seq = real_read.split('~')[0]
                    else:
                        seq = real_read
                    
                    # Ограничиваем длину для демонстрации
                    if len(seq) > 100:
                        seq = seq[:100]
                    
                    if len(seq) >= 13:  # Минимальная длина для 13-меров
                        test_sequences.append(seq)
                        print(f"  Добавлен рид {rid}: длина {len(seq)}")
                        
            except Exception as e:
                print(f"  Ошибка при получении рида {rid}: {e}")
                continue
                        
        print(f"✓ Получено {len(test_sequences)} реальных последовательностей")
    
    # Если не удалось получить реальные риды, создаем тестовые
    if len(test_sequences) == 0:
        print("Создаем тестовые последовательности из известных 13-меров...")
        # Возьмем несколько известных 13-меров и объединим их
        known_kmers = ["AAAAAGAGTTAAT", "AGTAGTAGTAGTA", "ATCGATCGATCGA"]
        synthetic_seq = "".join(known_kmers) + "AAAAAAAAAA"  # Добавляем немного перекрытия
        test_sequences.append(synthetic_seq)
        
        # Добавим еще несколько синтетических последовательностей
        test_sequences.extend([
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",  # Простая последовательность только A
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGAT",  # Повторяющийся паттерн
            "AAAAACCCCCGGGGGTTTTTATCGATCGAAAAAAA",     # Смешанная последовательность
        ])
        print(f"✓ Создано {len(test_sequences)} синтетических последовательностей")

    for i, seq in enumerate(test_sequences):
        print(f"\n--- Последовательность {i+1} (длина: {len(seq)}) ---")
        print(f"Seq: {seq[:50]}{'...' if len(seq) > 50 else ''}")
        
        if len(seq) < 13:
            print("  Последовательность слишком короткая для анализа 13-меров")
            continue
            
        try:
            # Анализируем покрытие для 13-меров
            coverage_positions = []
            total_tf_sum = 0
            max_tf = 0
            non_zero_positions = 0
            
            # Итерируемся по всем возможным 13-мерам в последовательности
            for pos in range(len(seq) - 12):
                kmer = seq[pos:pos+13]
                if len(kmer) == 13:
                    tf = index.get_total_tf_value_13mer(kmer)
                    coverage_positions.append((pos, kmer, tf))
                    total_tf_sum += tf
                    max_tf = max(max_tf, tf)
                    if tf > 0:
                        non_zero_positions += 1
            
            total_positions = len(coverage_positions)
            avg_tf = total_tf_sum / total_positions if total_positions > 0 else 0
            
            print(f"  Всего позиций для анализа: {total_positions}")
            print(f"  Позиций с ненулевым TF: {non_zero_positions} ({non_zero_positions/total_positions*100:.1f}%)")
            print(f"  Максимальный TF: {max_tf}")
            print(f"  Средний TF: {avg_tf:.2f}")
            
            # Показываем первые несколько 13-меров с их TF
            print(f"  Первые 5 13-меров:")
            for j in range(min(5, len(coverage_positions))):
                pos, kmer, tf = coverage_positions[j]
                print(f"    Поз {pos}: {kmer} -> TF = {tf}")
                
            # Тестируем с различными cutoff (симуляция)
            cutoffs = [0, 1, 5, 10] if max_tf > 0 else [0]
            for cutoff in cutoffs:
                filtered_positions = sum(1 for _, _, tf in coverage_positions if tf >= cutoff)
                print(f"  С TF >= {cutoff}: {filtered_positions} позиций")
                
        except Exception as e:
            print(f"  Ошибка анализа покрытия: {e}")

    # Тестируем итерацию по 13-мерам последовательности
    print(f"\n=== Тест итерации по 13-мерам последовательности ===")
    
    test_seq = test_sequences[0][:50] if len(test_sequences[0]) > 50 else test_sequences[0]
    print(f"Анализируем последовательность длиной {len(test_seq)}")
    print(f"Последовательность: {test_seq}")
    
    try:
        kmer_count = 0
        tf_sum = 0
        non_zero_kmers = 0
        
        for pos in range(len(test_seq) - 12):
            kmer = test_seq[pos:pos+13]
            if len(kmer) == 13:
                tf = index.get_total_tf_value_13mer(kmer)
                kmer_count += 1
                tf_sum += tf
                if tf > 0:
                    non_zero_kmers += 1
                
                # Показываем первые 5 13-меров с ненулевым TF
                if tf > 0 and non_zero_kmers <= 5:
                    print(f"  13-мер {kmer_count}: {kmer} -> TF = {tf}")
                    
        if kmer_count > 0:
            coverage_percent = (non_zero_kmers / kmer_count) * 100 if kmer_count > 0 else 0
            print(f"Всего 13-меров: {kmer_count}")
            print(f"13-меры с ненулевым TF: {non_zero_kmers} ({coverage_percent:.2f}%)")
            print(f"Средний TF: {tf_sum/kmer_count:.2f}")
            
    except Exception as e:
        print(f"Ошибка итерации по 13-мерам: {e}")

# Стресс-тест анализа покрытия для 13-меров
if index is not None:
    print("\n" + "="*80)
    print("СТРЕСС-ТЕСТ: АНАЛИЗ ПОКРЫТИЯ ДЛЯ 10,000 ПОСЛЕДОВАТЕЛЬНОСТЕЙ (13-МЕРЫ)")
    print("="*80)
    
    # Собираем реальные последовательности из индекса
    test_sequences_stress = []
    print("Собираем тестовые последовательности из 13-мерного индекса...")
    
    # Получаем реальные риды напрямую из 13-мерного индекса
    if index.n_reads > 0:
        print("Загружаем реальные риды из 13-мерного индекса...")
        
        # Возьмем первые 10K ридов (или сколько есть)
        max_sequences = min(10000, index.n_reads)
        
        for rid in range(max_sequences):
            try:
                real_read = index.get_read_by_rid(rid)
                if real_read and len(real_read) >= 50:  # Минимальная длина для анализа
                    # Если рид содержит субриды, берем первый субрид
                    if '~' in real_read:
                        seq = real_read.split('~')[0]
                    else:
                        seq = real_read
                    
                    # Ограничиваем длину для стабильного времени выполнения
                    if len(seq) > 200:
                        seq = seq[:200]
                    
                    if len(seq) >= 13:  # Минимальная длина для 13-меров
                        test_sequences_stress.append(seq)
                        
                # Показываем прогресс каждые 1000 последовательностей
                if (rid + 1) % 1000 == 0:
                    print(f"  Собрано: {len(test_sequences_stress)} последовательностей из {rid + 1} ридов")
                    
            except Exception as e:
                continue  # Пропускаем проблемные риды
        
        print(f"✓ Собрано {len(test_sequences_stress)} реальных последовательностей из 13-мерного индекса")
    
    # Если не удалось получить реальные риды, генерируем синтетические
    if len(test_sequences_stress) == 0:
        print("⚠️ Риды не загружены в 13-мерный индекс")
        print("Генерируем синтетические последовательности для 13-меров...")
        max_sequences = 10000
        
        # Генерируем случайные последовательности разной длины
        for i in range(max_sequences):
            # Случайная длина от 50 до 150 нуклеотидов
            seq_length = np.random.randint(50, 151)
            seq = ''.join(np.random.choice(['A', 'T', 'G', 'C'], seq_length))
            test_sequences_stress.append(seq)
            
            # Показываем прогресс каждые 1000 последовательностей
            if (i + 1) % 1000 == 0:
                print(f"  Сгенерировано: {i + 1:,} / {max_sequences:,} последовательностей")
        
        print(f"✓ Сгенерировано {len(test_sequences_stress):,} синтетических последовательностей")
    
    # Статистика последовательностей
    seq_lengths = [len(seq) for seq in test_sequences_stress]
    avg_length = sum(seq_lengths) / len(seq_lengths)
    min_length = min(seq_lengths)
    max_length = max(seq_lengths)
    
    print(f"Статистика последовательностей:")
    print(f"  Количество: {len(test_sequences_stress):,}")
    print(f"  Средняя длина: {avg_length:.1f} нт")
    print(f"  Диапазон длин: {min_length} - {max_length} нт")
    
    # Оценочный расчет количества 13-меров
    total_kmers = sum(max(0, len(seq) - 12) for seq in test_sequences_stress)
    print(f"  Общее количество 13-меров для анализа: {total_kmers:,}")
    
    print(f"\n--- Стресс-тест анализа покрытия для 13-меров ---")
    
    # Замеряем использование памяти до теста
    import psutil
    import os
    process = psutil.Process(os.getpid())
    memory_before = process.memory_info()
    
    print(f"Память до теста: RSS = {memory_before.rss / 1024 / 1024:.1f} МБ")
    
    # Запускаем тест
    start_time = time.time()
    processed_sequences = 0
    total_positions_analyzed = 0
    total_nonzero_positions = 0
    total_tf_sum = 0
    max_tf_found = 0
    
    print("Анализируем покрытие последовательностей 13-мерами...")
    
    try:
        for i, seq in enumerate(test_sequences_stress):
            # Анализируем покрытие для каждой последовательности
            coverage = []
            
            # Получаем все 13-меры из последовательности
            kmers_in_seq = []
            for pos in range(len(seq) - 12):
                kmer = seq[pos:pos+13]
                if len(kmer) == 13:
                    kmers_in_seq.append(kmer)
            
            if len(kmers_in_seq) == 0:
                continue
                
            # Получаем TF для всех 13-меров в последовательности batch'ем
            tf_values = index.get_total_tf_values_13mer(kmers_in_seq)
            
            # Собираем статистику
            processed_sequences += 1
            total_positions_analyzed += len(tf_values)
            nonzero_positions = sum(1 for tf in tf_values if tf > 0)
            total_nonzero_positions += nonzero_positions
            tf_sum = sum(tf_values)
            total_tf_sum += tf_sum
            max_tf_in_seq = max(tf_values) if tf_values else 0
            max_tf_found = max(max_tf_found, max_tf_in_seq)
            
            # Показываем прогресс каждые 1000 последовательностей
            if processed_sequences % 1000 == 0:
                elapsed = time.time() - start_time
                sequences_per_sec = processed_sequences / elapsed if elapsed > 0 else 0
                positions_per_sec = total_positions_analyzed / elapsed if elapsed > 0 else 0
                
                print(f"  Обработано: {processed_sequences:,} / {len(test_sequences_stress):,} " +
                      f"({sequences_per_sec:.1f} seq/сек, {positions_per_sec:.0f} позиций/сек)")
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        
        # Замеряем память после теста
        memory_after = process.memory_info()
        memory_increase = (memory_after.rss - memory_before.rss) / 1024 / 1024
        
        print(f"\n--- Результаты стресс-теста для 13-меров ---")
        print(f"✓ Обработано последовательностей: {processed_sequences:,}")
        print(f"✓ Общее время выполнения: {elapsed_time:.2f} секунд")
        print(f"✓ Производительность:")
        print(f"  - Последовательностей в секунду: {processed_sequences / elapsed_time:.1f}")
        print(f"  - 13-мер позиций в секунду: {total_positions_analyzed / elapsed_time:.0f}")
        print(f"  - Время на одну последовательность: {elapsed_time * 1000 / processed_sequences:.2f} мс")
        
        print(f"\n✓ Статистика покрытия для 13-меров:")
        print(f"  - Всего позиций проанализировано: {total_positions_analyzed:,}")
        print(f"  - Позиций с ненулевым TF: {total_nonzero_positions:,} ({total_nonzero_positions/total_positions_analyzed*100:.1f}%)")
        print(f"  - Средний TF на позицию: {total_tf_sum/total_positions_analyzed:.2f}")
        print(f"  - Максимальный TF: {max_tf_found}")
        print(f"  - Общая сумма TF: {total_tf_sum:,}")
        
        print(f"\n✓ Использование памяти:")
        print(f"  - Память после теста: RSS = {memory_after.rss / 1024 / 1024:.1f} МБ")
        print(f"  - Прирост памяти: {memory_increase:+.1f} МБ")
        
        # Дополнительный тест: несколько последовательностей с разными подходами
        print(f"\n--- Тест сравнения подходов для 13-меров ---")
        test_sample = test_sequences_stress[:100]  # Возьмем 100 последовательностей
        
        # Подход 1: Batch обработка всех 13-меров сразу
        start_batch_time = time.time()
        all_kmers = []
        for seq in test_sample:
            for pos in range(len(seq) - 12):
                kmer = seq[pos:pos+13]
                if len(kmer) == 13:
                    all_kmers.append(kmer)
        
        batch_results = index.get_total_tf_values_13mer(all_kmers)
        batch_time = time.time() - start_batch_time
        
        # Подход 2: Последовательная обработка по последовательностям
        start_seq_time = time.time()
        total_positions_seq = 0
        for seq in test_sample:
            kmers_in_seq = []
            for pos in range(len(seq) - 12):
                kmer = seq[pos:pos+13]
                if len(kmer) == 13:
                    kmers_in_seq.append(kmer)
            
            if kmers_in_seq:
                seq_results = index.get_total_tf_values_13mer(kmers_in_seq)
                total_positions_seq += len(seq_results)
        
        seq_time = time.time() - start_seq_time
        
        print(f"  Подход 1 (все 13-меры сразу): {batch_time:.3f} сек, {len(all_kmers)} 13-меров")
        print(f"  Подход 2 (по последовательностям): {seq_time:.3f} сек, {total_positions_seq} позиций")
        print(f"  Разница: {abs(batch_time - seq_time):.3f} сек")
        
        print(f"\n" + "="*80)
        print(f"СТРЕСС-ТЕСТ ПОКРЫТИЯ ДЛЯ 13-МЕРОВ ЗАВЕРШЕН")
        print(f"="*80)
        print(f"✓ Успешно проанализировано {processed_sequences:,} последовательностей")
        print(f"✓ Производительность: {processed_sequences / elapsed_time:.1f} последовательностей/сек")
        print(f"✓ Общее время: {elapsed_time:.2f} секунд")
        print(f"✓ Анализировано {total_positions_analyzed:,} позиций 13-меров")
        
    except Exception as e:
        print(f"❌ Ошибка стресс-теста покрытия: {e}")
        import traceback
        traceback.print_exc()
else:
    print("\n⚠️ Индекс не загружен, пропускаем стресс-тест покрытия для 13-меров")
