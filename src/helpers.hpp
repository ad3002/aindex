#ifndef HELP_FILE_H
#define HELP_FILE_H

#include <cstdint>
#include <cstddef>

void printProgressBar(double progress);
void printDoubleProgressBars(double progress1, double progress2);
void printTripleProgressBars(double progress1, double progress2, double progress3);
void printTripleProgressBars(double progress1, double progress2, double progress3, size_t a, size_t b, size_t c, size_t d, size_t e);

#endif