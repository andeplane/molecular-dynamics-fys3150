//#include "math/morton.h"

//uint64_t reinterpretDoubleAsUInt(double d) {
//  int const doubleSize = sizeof(double);
//  char* array = reinterpret_cast<char*>(&d);
//  uint64_t result = 0;
//  for (int i = 0; i < doubleSize; ++i) {
//    result += (uint64_t)array[i] << (8*i);
//  }

//  return result;
//}

//bool lessThanZOrderDouble(double *a, double *b) {
//  uint64_t j = 0;
//  uint64_t x = 0;

//  for (int i = 0; i < 3; ++i) {
//    uint64_t y = reinterpretDoubleAsUInt(a[i]) ^ reinterpretDoubleAsUInt(b[i]);
//    if (x < y && x < (x ^ y)) {
//      j = i;
//      x = y;
//    }
//  }
//  return (a[j] - b[j]) > 0;
//}
