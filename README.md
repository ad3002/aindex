# aindex

Perfect hash based index for genomic data.

## Installation

With pip:

```bash
pip install aindex2
```

If you want to install the package from source or you don't have pip version for your system, you can do so by running the following commands:

```bash
git clone https://github.com/ad3002/aindex.git
cd aindex
make
pip install .
```

This will create the necessary executables in the `bin` directory.

To uninstall:

```bash
pip uninstall aindex2
pip uninstall clean
```

To clean up the compiled files, run:

```
make clean
```

### Mac Compilation Command

Currently unsupported in Makefile. But you can try to compile the Python wrapper on MacOs manually with the following command:

```
g++ -c -std=c++11 -fPIC python_wrapper.cpp -o python_wrapper.o && g++ -c -std=c++11 -fPIC kmers.cpp kmers.hpp debrujin.cpp debrujin.hpp hash.cpp hash.hpp read.cpp read.hpp settings.hpp settings.cpp && g++ -shared -Wl,-install_name,python_wrapper.so -o python_wrapper.so python_wrapper.o kmers.o debrujin.o hash.o read.o settings.o
```

## Usage

Compute all binary arrays:

```bash
FASTQ1=./tests/raw_reads.101bp.IS350bp25_1.fastq
FASTQ2=./tests/raw_reads.101bp.IS350bp25_2.fastq
OUTPUT_PREFIX=tests/raw_reads.101bp.IS350bp25

time python ~/Dropbox/workspace/aindex/scripts/compute_aindex.py -i $FASTQ1,$FASTQ2 -t fastq -o $OUTPUT_PREFIX --lu 2 -P 30
```

## Usage from Python

You can simply run **demo.py** or:

```python
import aindex

prefix_path = "tests/raw_reads.101bp.IS350bp25"
index = aindex.get_aindex(prefix_path)

kmer = "A"*23
rkmer = "T"*23
kid = kmer2tf.get_kid_by_kmer(kmer)
print(kmer2tf.get_kmer_info_by_kid(kid))
print(kmer2tf[kmer], kid, kmer2tf.get_kmer_by_kid(kid), len(kmer2tf.pos(kmer)), kmer2tf.get_strand(kmer), kmer2tf.get_strand(rkmer))

pos = kmer2tf.pos(kmer)[0]
print(pos)

print(kmer2tf.get_kid_by_kmer(kmer), kmer2tf.get_kid_by_kmer(rkmer))

print(kmer2tf.get_hash_size())

print(kmer2tf.get_read(0, 123, 0))

print(kmer2tf.get_read(0, 123, 1))


k = 23
for p in kmer2tf.pos("GCAGCTCAGCAGGACGGCCAACC"):
  print(kmer2tf.get_read(p, p+k))
  break


print(kmer2tf["GCAGCTCAGCAGGACGGCCAACC"])

sequence = kmer2tf.get_read(0, 1023, 0)

for kmer, tf in kmer2tf.iter_sequence_kmers(sequence):
  print(kmer, tf)
  break


k = 23
sequence = "TAAGTTATTATTTAGTTAATACTTTTAACAATATTATTAAGGTATTTAAAAAATACTATTATAGTATTTAACATAGTTAAATACCTTCCTTAATACTGTTAAATTATATTCAATCAATACATATATAATATTATTAAAATACTTGATAAGTATTATTTAGATATTAGACAAATACTAATTTTATATTGCTTTAATACTTAATAAATACTACTTATGTATTAAGTAAATATTACTGTAATACTAATAACAATATTATTACAATATGCTAGAATAATATTGCTAGTATCAATAATTACTAATATAGTATTAGGAAAATACCATAATAATATTTCTACATAATACTAAGTTAATACTATGTGTAGAATAATAAATAATCAGATTAAAAAAATTTTATTTATCTGAAACATATTTAATCAATTGAACTGATTATTTTCAGCAGTAATAATTACATATGTACATAGTACATATGTAAAATATCATTAATTTCTGTTATATATAATAGTATCTATTTTAGAGAGTATTAATTATTACTATAATTAAGCATTTATGCTTAATTATAAGCTTTTTATGAACAAAATTATAGACATTTTAGTTCTTATAATAAATAATAGATATTAAAGAAAATAAAAAAATAGAAATAAATATCATAACCCTTGATAACCCAGAAATTAATACTTAATCAAAAATGAAAATATTAATTAATAAAAGTGAATTGAATAAAATTTTGAAAAAAATGAATAACGTTATTATTTCCAATAACAAAATAAAACCACATCATTCATATTTTTTAATAGAGGCAAAAGAAAAAGAAATAAACTTTTATGCTAACAATGAATACTTTTCTGTCAAATGTAATTTAAATAAAAATATTGATATTCTTGAACAAGGCTCCTTAATTGTTAAAGGAAAAATTTTTAACGATCTTATTAATGGCATAAAAGAAGAGATTATTACTATTCAAGAAAAAGATCAAACACTTTTGGTTAAAACAAAAAAAACAAGTATTAATTTAAACACAATTAATGTGAATGAATTTCCAAGAATAAGGTTTAATGAAAAAAACGATTTAAGTGAATTTAATCAATTCAAAATAAATTATTCACTTTTAGTAAAAGGCATTAAAAAAATTTTTCACTCAGTTTCAAATAATCGTGAAATATCTTCTAAATTTAATGGAGTAAATTTCAATGGATCCAATGGAAAAGAAATATTTTTAGAAGCTTCTGACACTTATAAACTATCTGTTTTTGAGATAAAGCAAGAAACAGAACCATTTGATTTCATTTTGGAGAGTAATTTACTTAGTTTCATTAATTCTTTTAATCCTGAAGAAGATAAATCTATTGTTTTTTATTACAGAAAAGATAATAAAGATAGCTTTAGTACAGAAATGTTGATTTCAATGGATAACTTTATGATTAGTTACACATCGGTTAATGAAAAATTTCCAGAGGTAAACTACTTTTTTGAATTTGAACCTGAAACTAAAATAGTTGTTCAAAAAAATGAATTAAAAGATGCACTTCAAAGAATTCAAACTTTGGCTCAAAATGAAAGAACTTTTTTATGCGATATGCAAATTAACAGTTCTGAATTAAAAATAAGAGCTATTGTTAATAATATCGGAAATTCTCTTGAGGAAATTTCTTGTCTTAAATTTGAAGGTTATAAACTTAATATTTCTTTTAACCCAAGTTCTCTATTAGATCACATAGAGTCTTTTGAATCAAATGAAATAAATTTTGATTTCCAAGGAAATAGTAAGTATTTTTTGATAACCTCTAAAAGTGAACCTGAACTTAAGCAAATATTGGTTCCTTCAAGATAATGAATCTTTACGATCTTTTAGAACTACCAACTACAGCATCAATAAAAGAAATAAAAATTGCTTATAAAAGATTAGCAAAGCGTTATCACCCTGATGTAAATAAATTAGGTTCGCAAACTTTTGTTGAAATTAATAATGCTTATTCAATATTAAGTGATCCTAACCAAAAGGAAAAATATGATTCAATGCTGAAAGTTAATGATTTTCAAAATCGCATCAAAAATTTAGATATTAGTGTTAGATGACATGAAAATTTCATGGAAGAACTCGAACTTCGTAAGAACTGAGAATTTGATTTTTTTTCATCTGATGAAGATTTCTTTTATTCTCCATTTACAAAAA"

test_kmer = "TAAGTTATTATTTAGTTAATACT"
right_kmer = "AGTTAATACTTTTAACAATATTA"

print "Task 1. Get kmer frequency"
raw_input("\nReady?")
for i in xrange(len(sequence)-k+1):
    kmer = sequence[i:i+k]
    print "Position %s kmer %s freq = %s" % (i, kmer, index[kmer])

print "Task 2. Iter read by read, print the first 20 reads"
raw_input("\nReady?")
for i, read in enumerate(index.iter_reads()):
    if i == 20:
        break
    print i, read

print "Task 3. Iter reads by kmer, returs (start, next_read_start, read, pos_if_uniq|None, all_poses)"
raw_input("\nReady?")
for read in iter_reads_by_kmer(test_kmer, index):
    print read

print "Task 4. Get distances in reads for two kmers, returns a list of (rid, left_kmer_pos, right_kmer_pos) tuples."
raw_input("\nReady?")
print get_left_right_distances(test_kmer, right_kmer, index)

print "Task 5. Get layout for kmer, returns (max_pos, reads, lefts, rights, rids, starts), for details see source code"
raw_input("\nReady?")
max_pos, reads, lefts, rights, rids, starts = get_layout_for_kmer(right_kmer, index)
print "Central layout:"
for read in reads:
    print read
print "Left flanks:"
print lefts
print "Right flanks:"
print rights

print "Task 6. Iter reads by sequence, returns (start, next_read_start, read, pos_if_uniq|None, all_poses)"
raw_input("\nReady?")
sequence = "AATATTATTAAGGTATTTAAAAAATACTATTATAGTATTTAACATA"
for read in iter_reads_by_sequence(sequence, index):
    print read
```