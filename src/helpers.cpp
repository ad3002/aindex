
#include <iostream>
#include <iomanip>

#include "helpers.hpp"

void printProgressBar(double progress) {
    int barWidth = 70;
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "%\r";
    if (progress == 1.0) {
        std::cout << std::endl;
    }
    std::cout.flush();
}

void printDoubleProgressBars(double progress1, double progress2) {
    int barWidth = 35; // Half of the previous bar width since we need space for two bars
    std::cout << "[";
    int pos1 = barWidth * progress1;
    int pos2 = barWidth * progress2;
    
    // First progress bar
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos1) std::cout << "=";
        else if (i == pos1) std::cout << ">";
        else std::cout << " ";
    }
    
    std::cout << "] ";
    std::cout << int(progress1 * 100.0) << "% [";

    // Second progress bar
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos2) std::cout << "=";
        else if (i == pos2) std::cout << ">";
        else std::cout << " ";
    }

    std::cout << "] " << int(progress2 * 100.0) << "%\r";
    
    // Move to the next line when both progress bars are complete
    if (progress1 == 1.0 && progress2 == 1.0) {
        std::cout << std::endl;
    }
    
    std::cout.flush();
}

// Function to print three progress bars
void printTripleProgressBars(double progress1, double progress2, double progress3) {
    int barWidth = 70; // Width for each progress bar
    std::cout << "\033[3A"; // Move cursor up 3 lines

    // First progress bar
    std::cout << "[";
    int pos1 = barWidth * progress1;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos1) std::cout << "=";
        else if (i == pos1) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress1 * 100.0) << "%\r";
    std::cout << std::endl; // Move to next line

    // Second progress bar
    std::cout << "[";
    int pos2 = barWidth * progress2;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos2) std::cout << "=";
        else if (i == pos2) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress2 * 100.0) << "%\r";
    std::cout << std::endl; // Move to next line

    // Third progress bar
    std::cout << "[";
    int pos3 = barWidth * progress3;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos3) std::cout << "=";
        else if (i == pos3) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress3 * 100.0) << "%\r";
    std::cout << std::endl; // Move to next line

    std::cout.flush();
}

void printTripleProgressBars(double progress1, double progress2, double progress3, size_t a, size_t b, size_t c, size_t d, size_t e) {
    int barWidth = 70; // Width for each progress bar
    std::cout << "\033[3A"; // Move cursor up 3 lines

    // First progress bar
    std::cout << "[";
    int pos1 = barWidth * progress1;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos1) std::cout << "=";
        else if (i == pos1) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress1 * 100.0) << "%" << " " << a << " " << e << " \r";
    std::cout << std::endl; // Move to next line

    // Second progress bar
    std::cout << "[";
    int pos2 = barWidth * progress2;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos2) std::cout << "=";
        else if (i == pos2) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress2 * 100.0) << "%" << " " << b << "\r";
    std::cout << std::endl; // Move to next line

    // Third progress bar
    std::cout << "[";
    int pos3 = barWidth * progress3;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos3) std::cout << "=";
        else if (i == pos3) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress3 * 100.0) << "%" << " " << c << "/" << d << "\r";
    std::cout << std::endl; // Move to next line

    std::cout.flush();
}