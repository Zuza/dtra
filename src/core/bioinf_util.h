#ifndef BIOINF_UTIL
#define BIOINF_UTIL

#include <cassert>
#include <ctime>
#include <cstdlib>

inline int baseToInt(char base, const int nValue = -1) {
  base = toupper(base);

  if (base == 'A') return 0;
  if (base == 'C') return 1;
  if (base == 'T') return 2;
  if (base == 'G') return 3;
  
  // by default, N will return random number [0..3]
  if (nValue == -1) return rand()%4;
  return nValue;
}

inline char intToBase(const int num) {
  if (num == 0) return 'A';
  if (num == 1) return 'C';
  if (num == 2) return 'T';
  if (num == 3) return 'G';
  
  // hopefully this doesn't happen
  assert(0 <= num && num < 4);
  return 0;
}

#endif
