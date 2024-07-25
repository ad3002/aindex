#ifndef HELP_FILE_H
#define HELP_FILE_H

#include <cstdint>
#include <cstddef>

void printProgressBar(double progress);
void printDoubleProgressBars(double progress1, double progress2);
void printTripleProgressBars(double progress1, double progress2, double progress3);
void printTripleProgressBars(double progress1, double progress2, double progress3, uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t e);

#endif