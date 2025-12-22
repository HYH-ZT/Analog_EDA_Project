#pragma once
#include <cstdlib>
#include <iostream>

#define OOB_ABORT(msg) do { \
  std::cerr << "\n[OOB_ABORT] " << msg << "\n"; \
  std::abort(); \
} while(0)

#define CHECK_IDX(idx, n, where) do { \
  if ((idx) < 0 || (idx) >= (n)) { \
    std::cerr << "\n[OOB] " << where \
              << " idx=" << (idx) << " size=" << (n) << "\n"; \
    std::abort(); \
  } \
} while(0)

#ifndef M_PI
constexpr double M_PI = 3.141592653589793238462643383279502884;
#endif
