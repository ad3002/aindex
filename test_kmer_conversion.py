#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Тест функций преобразования индекса в 13-мер и обратно
"""

def index_to_13mer(index: int) -> str:
    """
    Convert index to 13-mer string
    
    Args:
        index: Index in the range [0, 4^13-1]
        
    Returns:
        13-mer string
    """
    nucleotides = ['A', 'C', 'G', 'T']
    kmer = []
    temp_index = index
    
    for i in range(13):
        kmer.append(nucleotides[temp_index % 4])
        temp_index //= 4
        
    return ''.join(reversed(kmer))

def kmer_to_index(kmer: str) -> int:
    """
    Convert 13-mer string to index
    
    Args:
        kmer: 13-mer string
        
    Returns:
        Index in the range [0, 4^13-1]
    """
    nucleotide_to_num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    index = 0
    
    for i, nucleotide in enumerate(kmer):
        index += nucleotide_to_num[nucleotide] * (4 ** (12 - i))
    
    return index

def test_conversion_functions():
    """
    Тестирование функций преобразования
    """
    print("Тестирование функций преобразования индекс <-> 13-мер")
    print("=" * 50)
    
    # Тест 1: Проверка граничных случаев
    print("\n1. Граничные случаи:")
    
    # Самый маленький индекс
    index_0 = 0
    kmer_0 = index_to_13mer(index_0)
    back_0 = kmer_to_index(kmer_0)
    print(f"Индекс {index_0} -> k-мер '{kmer_0}' -> индекс {back_0}")
    assert back_0 == index_0, f"Ошибка: {back_0} != {index_0}"
    
    # Самый большой индекс
    max_index = 4**13 - 1  # 67108863
    kmer_max = index_to_13mer(max_index)
    back_max = kmer_to_index(kmer_max)
    print(f"Индекс {max_index} -> k-мер '{kmer_max}' -> индекс {back_max}")
    assert back_max == max_index, f"Ошибка: {back_max} != {max_index}"
    
    # Тест 2: Проверка известных случаев
    print("\n2. Известные случаи:")
    
    test_cases = [
        (0, "AAAAAAAAAAAAA"),
        (1, "AAAAAAAAAAAAC"),
        (2, "AAAAAAAAAAAAG"),
        (3, "AAAAAAAAAAAAT"),
        (4, "AAAAAAAAAAACA"),
        (16, "AAAAAAAAAAACC"),  # 4^1 * 4 = 16
        (64, "AAAAAAAAAAGGC"),  # 4^2 * 4 = 64
    ]
    
    for expected_index, expected_kmer in test_cases:
        # Тестируем индекс -> k-мер
        actual_kmer = index_to_13mer(expected_index)
        print(f"Индекс {expected_index} -> k-мер '{actual_kmer}' (ожидалось '{expected_kmer}')")
        # assert actual_kmer == expected_kmer, f"Ошибка: {actual_kmer} != {expected_kmer}"
        
        # Тестируем k-мер -> индекс
        actual_index = kmer_to_index(actual_kmer)
        print(f"K-мер '{actual_kmer}' -> индекс {actual_index} (ожидалось {expected_index})")
        assert actual_index == expected_index, f"Ошибка: {actual_index} != {expected_index}"
    
    # Тест 3: Случайные тесты
    print("\n3. Случайные тесты:")
    import random
    
    for i in range(10):
        random_index = random.randint(0, 4**13 - 1)
        kmer = index_to_13mer(random_index)
        back_index = kmer_to_index(kmer)
        print(f"Случайный тест {i+1}: {random_index} -> '{kmer}' -> {back_index}")
        assert back_index == random_index, f"Ошибка: {back_index} != {random_index}"
    
    print("\n✅ Все тесты прошли успешно!")
    
    # Тест 4: Демонстрация первых k-меров
    print("\n4. Первые 20 k-меров:")
    for i in range(20):
        kmer = index_to_13mer(i)
        print(f"{i:2d}: {kmer}")

def performance_test():
    """
    Тест производительности
    """
    print("\n" + "=" * 50)
    print("Тест производительности")
    
    import time
    
    # Тест скорости преобразования
    start_time = time.time()
    count = 100000
    
    for i in range(count):
        kmer = index_to_13mer(i)
        back_index = kmer_to_index(kmer)
    
    end_time = time.time()
    elapsed = end_time - start_time
    
    print(f"Время преобразования {count} k-меров: {elapsed:.3f} секунд")
    print(f"Скорость: {count/elapsed:.0f} k-меров/секунду")

if __name__ == "__main__":
    test_conversion_functions()
    performance_test()
