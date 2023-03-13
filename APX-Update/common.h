#pragma once

using namespace std;

#include <stdio.h>
#include <iostream>
#include <list>
#include <map>
#include <utility>
#include <vector>
#include <memory.h>
#include <algorithm>
#include <chrono>
#include <execinfo.h>
#include <numeric>

#define ASSERT(truth) \
    if (!(truth)) { \
        void *buffer[100]; \
        int nptrs = backtrace(buffer, 100); \
        char **strings = backtrace_symbols(buffer, nptrs); \
        if (strings) { \
            for (int i = 0; i < nptrs; ++i) { \
                printf("%s\n", strings[i]); \
            } \
            free(strings); \
        } \
        printf("ASSERT LINE %d, %s\n", __LINE__, __FILE__);\
      std::exit(EXIT_FAILURE); \
    } else

#define ASSERT_MSG(truth, msg) \
    if (!(truth)) { \
      std::cerr << "\x1b[1;31mASSERT\x1b[0m: " \
                << "LINE " << __LINE__ \
                << ", " << __FILE__ << '\n' \
                << "\x1b[1;32mINFO\x1b[0m: " << msg \
                << std::endl; \
      std::exit(EXIT_FAILURE); \
    } else


