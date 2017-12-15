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
```






