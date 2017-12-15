# aindex
Perfect hash based Index for text data

## Installation

```
cd external
git clone https://github.com/ot/emphf.git
cd emphf
cmake .
make
```

```
g++ Compute_index.cpp  -std=c++11 -pthread -O3 debrujin.hpp debrujin.cpp read.hpp read.cpp kmers.hpp kmers.cpp settings.cpp settings.hpp hash.hpp hash.cpp && mv a.out bin/compute_index.exe

g++ Compute_aindex.cpp  -std=c++11 -pthread -O3 debrujin.hpp debrujin.cpp read.hpp read.cpp kmers.hpp kmers.cpp settings.cpp settings.hpp hash.hpp hash.cpp && mv a.out bin/compute_aindex.exe

g++ Compute_reads.cpp  -std=c++11 -pthread -O3 debrujin.hpp debrujin.cpp read.hpp read.cpp kmers.hpp kmers.cpp settings.cpp settings.hpp hash.hpp hash.cpp && mv a.out bin/compute_reads.exe
```

Lunux compilation commadn:

```
g++ -c -std=c++11 -fPIC python_wrapper.cpp -o python_wrapper.o && g++ -c -std=c++11 -fPIC kmers.cpp kmers.hpp debrujin.cpp debrujin.hpp hash.cpp hash.hpp read.cpp read.hpp settings.hpp settings.cpp && g++ -shared -Wl,-soname,python_wrapper.so -o python_wrapper.so python_wrapper.o kmers.o debrujin.o hash.o read.o settings.o
```

Mac compilation command:

```
g++ -c -std=c++11 -fPIC python_wrapper.cpp -o python_wrapper.o && g++ -c -std=c++11 -fPIC kmers.cpp kmers.hpp debrujin.cpp debrujin.hpp hash.cpp hash.hpp read.cpp read.hpp settings.hpp settings.cpp && g++ -shared -Wl,-install_name,python_wrapper.so -o python_wrapper.so python_wrapper.o kmers.o debrujin.o hash.o read.o settings.o
```

## Usage


Compute all binary arrays:

```bash
FASTQ1=raw_reads.101bp.IS350bp25_1.fastq
FASTQ2=raw_reads.101bp.IS350bp25_2.fastq

jellyfish count -m 23 -s 2G -t 4 -C -o kmers.23.jf2 $FASTQ1 $FASTQ2
jellyfish dump -c -t -o kmers.23.dat kmers.23.jf2
cut -f1 kmers.23.dat > kmer.23.kmers
../external/emphf/compute_mphf_seq kmer.23.kmers kmer.23.p

../bin/compute_index.exe kmers.23.dat kmers.23.pf kemrs.23 4 0

../bin/compute_reads.exe raw_reads.101bp.IS350bp25_1.fastq raw_reads.101bp.IS350bp25_2.fastq fastq reads.reads

../bin/compute_aindex.exe reads.reads kmers.23.pf kmers.23 kmers.23 4 23 kmers.23.tf.bin 

```











