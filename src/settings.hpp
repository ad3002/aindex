//
// Created by Aleksey Komissarov on 02/09/15.
//

#ifndef STIRKA_SETTINGS_H
#define STIRKA_SETTINGS_H

#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <locale>
#include <memory>
#include <ctime>
#include <cstring>

namespace Settings {

    extern uint32_t MINIMAL_READ_LENGTH;
    extern uint32_t TRUE_ERRORS;
    extern uint32_t TRUE_NOT_ERRORS;
    extern uint32_t TRUE_REF_ERRORS;
    extern uint32_t PRE_ERRORS;
    extern uint32_t K;
    extern uint32_t JFK;
    extern uint32_t MIN_Q;
    extern uint32_t MIN_MI;
    extern uint32_t MIN_MI_N;
    extern uint32_t MAX_MA;
    extern uint32_t TRIM_LOW_COV_TOLERANCE;
    extern uint32_t STIRKA_FILL_OVERLAP;
    extern uint32_t STIRKA_RESOLVE_SNP;
    extern uint32_t VERBOSE;
    extern uint32_t CORRECTION_LEVEL;
    extern uint32_t COVERAGE;
    extern uint32_t MERGE;
    extern uint32_t INDEX_SIZE;

}
#endif //STIRKA_SETTINGS_H
