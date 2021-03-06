#ifndef BIOINF_UTIL
#define BIOINF_UTIL

#include <cassert>
#include <cctype>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <random>

inline bool isBase(char isit) {
  isit = toupper(isit);

  return isit == 'A' || isit == 'T' || isit == 'G' || isit == 'C';
}


// translate bases (ie. 'Y' -> C or T), ref: http://www.bioinformatics.org/sms/iupac.html
inline char randBaseToBase(char base) {
  thread_local int len[26] = {1, 3, 1, 3, 0, 0, 1, 3, 0, 0, 2, 0, 2, 4, 0, 0, 0, 2, 2, 1, 1, 3, 2, 0, 2, 0};
  thread_local char trans[26][5] = {"A",
                              "CGT",
                              "C",
                              "AGT",
                              "",
                              "",
                              "G",
                              "ACT",
                              "",
                              "",
                              "GT",
                              "",
                              "AC",
                              "AGTC",
                              "",
                              "",
                              "",
                              "AG",
                              "GC",
                              "T",
                              "U",
                              "ACG",
                              "AT",
                              "",
                              "CT",
                              ""};

  assert(base >= 'A' && base <= 'Z');

  if (len[base-'A'] == 1) {
    return base;
  }
  
  int t = len[base-'A'];
  assert(t != 0);
  t = rand() % t;
  return trans[base-'A'][t];
}


inline int baseToInt(char base, const int nValue = -1) {
  base = toupper(base);

  if (base == 'A') return 0;
  if (base == 'C') return 1;
  if (base == 'T') return 2;
  if (base == 'G') return 3;

  base = randBaseToBase(base);

  if (base == 'A') return 0;
  if (base == 'C') return 1;
  if (base == 'T') return 2;
  if (base == 'G') return 3;

  assert(false);
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

// substitute every occurence of 'N' with another random
inline void subNonAcgtWithRandom(char* data, size_t data_len) {
  for (size_t i = 0; i < data_len; ++i) {
    if (data[i] != 'A' && data[i] != 'C' && data[i] != 'G' && data[i] != 'T') {
      data[i] = randBaseToBase(data[i]);
    }
  }
}

inline char getBaseComplement(char base) {
  base = toupper(base);

  if (base == 'A') return 'T';
  if (base == 'T') return 'A';
  if (base == 'C') return 'G';
  if (base == 'G') return 'C';

  if (base != 'N') {
    fprintf(stderr, "invalid base in getBaseComplement: %c\n", base);
    return base;
  }
  assert(base == 'N');
  return 'N';
}


inline std::string ReverseComplement(std::string s) {
  std::string ret;
  ret.reserve(s.size()+1);
  for (int i = (int)s.size()-1; i >= 0; --i) {
    if (s[i] == 'A') ret += "T";
    if (s[i] == 'T') ret += "A";
    if (s[i] == 'C') ret += "G";
    if (s[i] == 'G') ret += "C";
    if (s[i] == 'N') ret += "N";
  }
  return ret;
}

#endif
