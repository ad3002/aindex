# aindex
Perfect hash based Index for text data

## Installation

```
git clone https://github.com/ot/emphf.git
cd emphf
cmake .
make
```

```
g++ Compute_index.cpp  -std=c++11 -pthread -O3 debrujin.hpp debrujin.cpp read.hpp read.cpp kmers.hpp kmers.cpp settings.cpp settings.hpp hash.hpp hash.cpp 
mv a.out bin/compute_index.exe

g++ Compute_aindex.cpp  -std=c++11 -pthread -O3 debrujin.hpp debrujin.cpp read.hpp read.cpp kmers.hpp kmers.cpp settings.cpp settings.hpp hash.hpp hash.cpp 
mv a.out bin/compute_aindex.exe
```

Lunux compilation commadn:

```
g++ -c -std=c++11 -fPIC python_wrapper.cpp -o python_wrapper.o && g++ -c -std=c++11 -fPIC kmers.cpp kmers.hpp debrujin.cpp debrujin.hpp hash.cpp hash.hpp read.cpp read.hpp settings.hpp settings.cpp && g++ -shared -Wl,-soname,python_wrapper.so -o python_wrapper.so python_wrapper.o kmers.o debrujin.o hash.o read.o settings.o
```

Mac compilation command:

```
g++ -c -std=c++11 -fPIC python_wrapper.cpp -o python_wrapper.o && g++ -c -std=c++11 -fPIC kmers.cpp kmers.hpp debrujin.cpp debrujin.hpp hash.cpp hash.hpp read.cpp read.hpp settings.hpp settings.cpp && g++ -shared -Wl,-install_name,python_wrapper.so -o python_wrapper.so python_wrapper.o kmers.o debrujin.o hash.o read.o settings.o
```











